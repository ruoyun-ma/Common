package kernel;

import rs2d.commons.log.Log;
import rs2d.spinlab.instrument.Instrument;
import rs2d.spinlab.instrument.util.GradientMath;
import rs2d.spinlab.sequence.SequenceTool;
import rs2d.spinlab.sequence.element.Opcode;
import rs2d.spinlab.sequence.table.Table;
import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.sequenceGenerator.util.TimeEvents;
import rs2d.spinlab.tools.param.*;
import rs2d.spinlab.tools.table.Order;

import java.util.*;

import static common.CommonSP.Grad_shape_rise_down;
import static common.CommonSP.Grad_shape_rise_up;
import static common.CommonSP.Time_grad_ramp;
import static java.util.Arrays.asList;

import common.*;
import model.*;

import static common.CommonUP.*;
import static common.CommonSP.*;
import static kernel.KernelGE.UP.*;
import static kernel.KernelGE.SP.*;

/**
 * Abstract Class KernelGE
 * prep Kernel functions of Gradient Echo Sequences
 * V1.0- 2021-4-26 XG
 * <p>
 * Keyhole/ TOF/ Flow/ Elliptical Sampling etc. are kept in each individual sequence
 * <p>
 */

public abstract class KernelGE extends SeqPrep {
    protected String sequenceVersion = "Version KernelGE 1.0";
    protected boolean CameleonVersion105 = false;

    protected String kspace_filling;
    protected double min_flush_delay;

    protected double slice_thickness_excitation;

    protected double txLength90;
    protected List<Double> grad_amp_spoiler_sl_ph_re;

    //Dynamic
    protected boolean isDynamic;
    protected int nb_dynamic_acquisition;
    protected boolean isDynamicMinTime;
    protected double time_between_frames;

    protected boolean is_rf_spoiling;

    protected boolean is_grad_rewinding;
    protected boolean is_grad_spoiler;
    protected boolean is_k_s_centred;

    protected double grad_shape_rise_time;
    protected RFPulse pulseRX;
    protected Gradient gradPhase2D;
    protected Gradient gradReadPrep;
    protected Gradient gradReadout;
    protected Gradient gradSliceRefPhase3D;

    protected enum UP implements GeneratorParamEnum {
        KSPACE_FILLING,
        KS_CENTERED,
        DYNAMIC_SEQUENCE,
        DYN_NUMBER_OF_ACQUISITION,
        DYNAMIC_MIN_TIME,
        DYN_TIME_BTW_FRAMES,
        RF_SPOILING,
        GRADIENT_ENABLE_REWINDING,
        GRADIENT_ENABLE_SPOILER,
        GRAD_AMP_SPOILER_SL_PH_RE,
        SLICE_REFOCUSING_GRADIENT_RATIO,
        NUMBER_OF_INTERLEAVED_SLICE,
        GRADIENT_PHASE_APPLICATION_TIME,
        GRADIENT_SPOILER_TIME,
        TRIGGER_TIME,

        //Below UP maybe removed in next version
        INTERLEAVED_ECHO_TRAIN,
        INTERLEAVED_NUM_OF_ECHO_TRAIN,
        ;

        @Override
        public Param build() {
            return null;
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum {
        Time_rf_spoiling,
        FreqOffset_RFSpoiling,

        Grad_enable_read,
        Grad_enable_phase_2D,
        Grad_enable_phase_3D,
        Grad_enable_slice,
        Grad_enable_spoiler_slice,
        Grad_enable_spoiler_phase,
        Grad_enable_spoiler_read,
        Grad_amp_read_read,
        FreqOffset_tx_prep,
        FreqOffset_tx_comp,
        Frequency_offset_init,
        Time_btw_dyn_frames,

        Grad_amp_spoiler_slice,
        Grad_amp_spoiler_phase,
        Grad_amp_spoiler_read,
        Time_grad_spoiler_top,

        //Below UP maybe removed in next version
        ;
    }

    public KernelGE() {
        super();
        addUserParams();
    }

    @Override
    public void init() {
        super.init();
        ((TextParam) getParam(TX_SHAPE)).setSuggestedValues(tx_shape);
        ((TextParam) getParam(TX_SHAPE)).setRestrictedToSuggested(true);

        //TRANSFORM PLUGIN
        TextParam transformPlugin = getParam(CommonUP.TRANSFORM_PLUGIN);
        transformPlugin.setSuggestedValues(asList("Sequential4D", "Sequential4DBackAndForth", "EPISequential4D", "Centric4D"));
        transformPlugin.setRestrictedToSuggested(true);


        // KSPACE_FILLING
        TextParam ksFilling = getParam(KSPACE_FILLING);
        ksFilling.setSuggestedValues(asList("Linear"));
        ksFilling.setRestrictedToSuggested(true);
    }

    // ==============================
// -----   GENERATE
// ==============================
    @Override
    public void initUserParam() {
        super.initUserParam();
        getParam(SEQUENCE_VERSION).setValue(sequenceVersion);

        kspace_filling = getText(KSPACE_FILLING);
        txLength90 = getDouble(TX_LENGTH);

        isDynamic = getBoolean(DYNAMIC_SEQUENCE);
        nb_dynamic_acquisition = isDynamic ? getInt(DYN_NUMBER_OF_ACQUISITION) : 1;
        isDynamic = isDynamic && (nb_dynamic_acquisition > 1);
        isDynamicMinTime = getBoolean(DYNAMIC_MIN_TIME);

        is_rf_spoiling = getBoolean(RF_SPOILING);

        is_grad_rewinding = getBoolean(GRADIENT_ENABLE_REWINDING);// get slice refocussing ratio
        is_grad_spoiler = getBoolean(GRADIENT_ENABLE_SPOILER);// get slice refocussing ratio
        grad_amp_spoiler_sl_ph_re = getListDouble(GRAD_AMP_SPOILER_SL_PH_RE);
        is_k_s_centred = getBoolean(KS_CENTERED);
    }

    //--------------------------------------------------------------------------------------
    // ini functions for beforeRouting() and initAfterRouting()
    //--------------------------------------------------------------------------------------
    @Override
    protected void iniModels() {
        setModels(new ArrayList<>(Arrays.asList("ExtTrig")), this);
    }

    @Override
    protected void iniTransformPlugin() throws Exception {
        ///// transform plugin
        if (isMultiplanar) {
            kspace_filling = "Linear";
            getParam(KSPACE_FILLING).setValue(kspace_filling);
            getParam(TRANSFORM_PLUGIN).setValue("Sequential4D");
        } else {
            switch (kspace_filling) {
                default:
                    kspace_filling = "Linear";
                    getParam(KSPACE_FILLING).setValue(kspace_filling);
                    getParam(TRANSFORM_PLUGIN).setValue("Sequential4D");
                    break;
            }
        }

        plugin = getTransformPlugin();
        plugin.setParameters(new ArrayList<>(getUserParams()));
        plugin.getScanOrder(); //output traj. graphs for Elliptical3D plugin
    }

    @Override
    protected void iniScanLoop() {
        nb_scan_1d = nb_averages;
        String updateDimension = SequenceTool.UpdateDimensionEnum.Disable.getText();
        if (!isMultiplanar) {
            nb_scan_2d = acqMatrixDimension2D;
            nb_scan_3d = acqMatrixDimension3D;
        } else {
            nb_scan_2d = acqMatrixDimension2D;
            nb_scan_3d = nb_shoot_3d;
        }

        //Dynamic and multi echo are filled into the 4th Dimension for 2D Imaging
        if (isMultiplanar) {
            nb_scan_4d = ((ExtTrig) models.get("ExtTrig")).nb_trigger * nb_dynamic_acquisition;
        } else {
            nb_scan_4d = Math.max(acqMatrixDimension4D / nb_interleaved_slice, 1);
        }

        nb_planar_excitation = (isMultiplanar ? acqMatrixDimension3D : 1);

        if (isKSCenterMode) { // Do only the center of the k-space for auto RG
            nb_scan_1d = 1;
            nb_scan_2d = 2;
            nb_scan_3d = !isMultiplanar ? 1 : nb_scan_3d;
            nb_scan_4d = 1;
            getParam(ACQUISITION_MATRIX_DIMENSION_3D).setValue(!isMultiplanar ? 1 : acqMatrixDimension3D);
            getParam(ACQUISITION_MATRIX_DIMENSION_4D).setValue(1);
        }
        // set Nb_scan  Values
        set(Pre_scan, nb_preScan); // Do the prescan
        set(Nb_point, acqMatrixDimension1D);
        set(Nb_1d, nb_scan_1d);
        set(Nb_2d, nb_scan_2d);
        set(Nb_3d, nb_scan_3d);
        //set(Update_dimension, updateDimension);
        set(Nb_4d, nb_scan_4d);
        set(Nb_interleaved_slice, nb_interleaved_slice - 1);
    }

    @Override
    protected void iniSeqDisp() {
        String seqDescription = "GE_";
        if (isMultiplanar) {
            seqDescription += "2D_";
        } else {
            seqDescription += "3D_";
        }
        String orientation = getText(ORIENTATION);
        seqDescription += orientation.substring(0, 3);

        String seqMatrixDescription = "_";
        seqMatrixDescription += userMatrixDimension1D + "x" + acqMatrixDimension2D + "x" + acqMatrixDimension3D;
        if (acqMatrixDimension4D != 1) {
            seqMatrixDescription += "x" + acqMatrixDimension4D;
        }
        seqDescription += seqMatrixDescription;

        if (isDynamic) {
            seqDescription += "_DYN=" + nb_dynamic_acquisition;
        }
        getParam(SEQ_DESCRIPTION).setValue(seqDescription);
    }

    @Override
    protected void iniSeqParamBasics() {
        set(Time_min_instruction, minInstructionDelay);
        set(Phase_reset, PHASE_RESET);
        set(Frequency_offset_init, 0.0);// PSD should start with a zero offset frequency pulse
    }

    @Override
    protected void iniSeqParamEnabled() {
        set(Grad_enable_read, isEnableRead);              // pass gradient line status to sequence
        set(Grad_enable_phase_2D, isEnablePhase);
        set(Grad_enable_phase_3D, (!isMultiplanar && isEnablePhase3D) || isEnableSlice);
        set(Grad_enable_slice, isEnableSlice);

        set(Grad_enable_spoiler_slice, (((!isMultiplanar && is_grad_rewinding && isEnablePhase3D) || (is_grad_rewinding && isEnableSlice) || (is_grad_spoiler && (grad_amp_spoiler_sl_ph_re.get(0) != 0)))));
        set(Grad_enable_spoiler_phase, (isEnablePhase && (is_grad_rewinding) || (is_grad_spoiler && (grad_amp_spoiler_sl_ph_re.get(1) != 0))));
        set(Grad_enable_spoiler_read, (isEnableRead && (is_grad_rewinding) || (is_grad_spoiler && (grad_amp_spoiler_sl_ph_re.get(2) != 0))));
    }

    //--------------------------------------------------------------------------------------
    // prep functions for afterRouting()
    //--------------------------------------------------------------------------------------
    @Override
    protected void prepRFandSliceGrad() throws Exception {
        getGradRiseTime();
        // -----------------------------------------------
        // Calculation RF pulse parameters  1/4 : Pulse declaration & Fatsat Flip angle calculation
        // -----------------------------------------------
        set(Time_tx, txLength90);    // set RF pulse length to sequence
        pulseTX = RFPulse.createRFPulse(getSequence(), Tx_att, Tx_amp, Tx_phase, Time_tx, Tx_shape, Tx_shape_phase, Tx_freq_offset);
        double flip_angle = getDouble(FLIP_ANGLE);

        // -----------------------------------------------
        // Calculation RF pulse parameters  2/4 : Shape
        // -----------------------------------------------
        pulseTX.setShape((getText(TX_SHAPE)), nb_shape_points, "Hamming");

        // -----------------------------------------------
        // Calculation RF pulse parameters  3/4 : RF pulse & attenuation
        // -----------------------------------------------
        if (getBoolean(TX_AMP_ATT_AUTO)) {
            if (!pulseTX.checkPower(flip_angle, observeFrequency, nucleus)) {
                notifyOutOfRangeParam(TX_LENGTH, pulseTX.getPulseDuration(), ((NumberParam) getParam(TX_LENGTH)).getMaxValue(), "Pulse length too short to reach RF power with this pulse shape");
                txLength90 = pulseTX.getPulseDuration();
            }
            pulseTX.prepAtt(80, getListInt(TX_ROUTE));
            pulseTX.prepTxAmp(getListInt(TX_ROUTE));
            rfPulses.add(pulseTX);
            rfPulsesTree.put(pulseTX.getPower(), pulseTX);
            getUPDisp();
        } else {
            pulseTX.setAtt(getInt(TX_ATT));
            pulseTX.setAmp(getDouble(TX_AMP));
            this.getParam(TX_AMP_90).setValue(getDouble(TX_AMP) * 90 / flip_angle);     // display 90° amplitude
            this.getParam(TX_AMP_180).setValue(getDouble(TX_AMP) * 90 / flip_angle);   // display 180° amplitude
        }

        // -----------------------------------------------
        // Calculation RF pulse parameters  4/4: bandwidth
        // -----------------------------------------------
        double tx_bandwidth_factor_90 = getTx_bandwidth_factor(TX_SHAPE, TX_BANDWIDTH_FACTOR, TX_BANDWIDTH_FACTOR_3D);
        double tx_bandwidth_90 = tx_bandwidth_factor_90 / txLength90;

        // ---------------------------------------------------------------------
        // calculate SLICE gradient amplitudes for RF pulses
        // ---------------------------------------------------------------------
        slice_thickness_excitation = (isMultiplanar ? sliceThickness : (sliceThickness * userMatrixDimension3D));
        gradSlice = Gradient.createGradient(getSequence(), Grad_amp_slice, Time_tx, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
        if (isEnableSlice && !gradSlice.prepareSliceSelection(tx_bandwidth_90, slice_thickness_excitation)) {
            slice_thickness_excitation = gradSlice.getSliceThickness();
            double slice_thickness_min = (isMultiplanar ? slice_thickness_excitation : (slice_thickness_excitation / userMatrixDimension3D));
            notifyOutOfRangeParam(SLICE_THICKNESS, slice_thickness_min, ((NumberParam) getParam(SLICE_THICKNESS)).getMaxValue(), "Pulse length too short to reach this slice thickness");
            sliceThickness = slice_thickness_min;
        }
        gradSlice.applyAmplitude();

        // ------------------------------------------------------------------
        //calculate TX FREQUENCY offsets tables for slice positionning
        // ------------------------------------------------------------------
        if (isMultiplanar && nb_planar_excitation > 1 && isEnableSlice) {
            //MULTI-PLANAR case : calculation of frequency offset table
            pulseTX.prepareOffsetFreqMultiSlice(gradSlice, nb_planar_excitation, spacingBetweenSlice, off_center_distance_3D);
            pulseTX.reoderOffsetFreq(plugin, acqMatrixDimension1D, nb_slices_acquired_in_single_scan);
            pulseTX.setFrequencyOffset(nb_slices_acquired_in_single_scan != 1 ? Order.ThreeLoop : Order.Three);
        } else {
            //3D CASE :
            pulseTX.prepareOffsetFreqMultiSlice(gradSlice, 1, 0, off_center_distance_3D);
            pulseTX.setFrequencyOffset(Order.Three);
        }

        // ------------------------------------------------------------------
        // calculate TX FREQUENCY offsets compensation
        // ------------------------------------------------------------------
        double grad_ratio_slice_refoc = isEnableSlice ? getDouble(SLICE_REFOCUSING_GRADIENT_RATIO) : 0.0;   // get slice refocussing ratio

        RFPulse pulseTXPrep = RFPulse.createRFPulse(getSequence(), Time_grad_ramp, FreqOffset_tx_prep);
        pulseTXPrep.setCompensationFrequencyOffset(pulseTX, grad_ratio_slice_refoc);

        RFPulse pulseTXComp = RFPulse.createRFPulse(getSequence(), Time_grad_ramp, FreqOffset_tx_comp);
        pulseTXComp.setCompensationFrequencyOffset(pulseTX, grad_ratio_slice_refoc);
    }

    @Override
    protected void prepDicom() {
        // Set  TRIGGER_TIME for dynamic or trigger acquisition
        if (isDynamic && (nb_dynamic_acquisition != 1) && !models.get("ExtTrig").isEnabled()) {
            ArrayList<Number> arrayListTrigger = new ArrayList<>();
            for (int i = 0; i < nb_dynamic_acquisition; i++) {
                arrayListTrigger.add(i * time_between_frames);
            }
            getParam(TRIGGER_TIME).setValue(arrayListTrigger);
        }
    }

    //--------------------------------------------------------------------------------------
    // get functions
    //--------------------------------------------------------------------------------------
    @Override
    protected void getAcq3D() {
        // MATRIX
        boolean is_partial_oversampling = getBoolean(PARTIAL_OVERSAMPLING);
        is_partial_oversampling = !isMultiplanar && userMatrixDimension3D >= 8 && is_partial_oversampling;
        getParam(PARTIAL_OVERSAMPLING).setValue(is_partial_oversampling);

        // 3D ZERO FILLING
        double partial_slice;
        if (isMultiplanar || is_partial_oversampling) {
            getParam(USER_PARTIAL_SLICE).setValue(100);
            partial_slice = 100;
        } else {
            partial_slice = getDouble(USER_PARTIAL_SLICE);
        }
        double zero_filling_3D = (100 - partial_slice) / 100f;
        getParam(USER_ZERO_FILLING_3D).setValue((100 - partial_slice));

        //Calculate the number of k-space lines acquired in the 3rd Dimension : acqMatrixDimension3D
        if (!isMultiplanar) {
            acqMatrixDimension3D = floorEven((1 - zero_filling_3D) * userMatrixDimension3D);
            acqMatrixDimension3D = (acqMatrixDimension3D < 4) && isEnablePhase3D ? 4 : acqMatrixDimension3D;
            userMatrixDimension3D = Math.max(userMatrixDimension3D, acqMatrixDimension3D);
            getParam(USER_MATRIX_DIMENSION_3D).setValue(userMatrixDimension3D);
        } else {
            if ((userMatrixDimension3D * 3 + ((is_rf_spoiling) ? 1 : 0) + 3 + 1) >= offset_channel_memory) {
                userMatrixDimension3D = ((int) Math.floor((offset_channel_memory - 4 - ((is_rf_spoiling) ? 1 : 0)) / 3.0));
                getParam(USER_MATRIX_DIMENSION_3D).setValue(userMatrixDimension3D);
            }
            acqMatrixDimension3D = userMatrixDimension3D;
        }

        if (!isMultiplanar) {
            nb_interleaved_slice = getInferiorDivisorToGetModulusZero(getInt(NUMBER_OF_INTERLEAVED_SLICE), acqMatrixDimension4D);
            nb_shoot_3d = acqMatrixDimension3D;
        } else {
            nb_shoot_3d = getInferiorDivisorToGetModulusZero(nb_shoot_3d, acqMatrixDimension3D);
            nb_interleaved_slice = (int) Math.ceil((acqMatrixDimension3D / (double) nb_shoot_3d));
        }

        getParam(NUMBER_OF_SHOOT_3D).setValue(nb_shoot_3d);
        getParam(NUMBER_OF_INTERLEAVED_SLICE).setValue(isMultiplanar ? nb_interleaved_slice : 0);

        acqMatrixDimension3D = is_partial_oversampling ? (int) Math.round(acqMatrixDimension3D / 0.8 / 2) * 2 : acqMatrixDimension3D;
        userMatrixDimension3D = is_partial_oversampling ? (int) Math.round(userMatrixDimension3D / 0.8 / 2) * 2 : userMatrixDimension3D;

        getParam(ACQUISITION_MATRIX_DIMENSION_3D).setValue(acqMatrixDimension3D);
        getParam(USER_MATRIX_DIMENSION_3D).setValue(userMatrixDimension3D);

        if (isMultiplanar) {
            spacingBetweenSlice = Math.max(getDouble(SPACING_BETWEEN_SLICE), 0);
            fov3d = sliceThickness * userMatrixDimension3D + spacingBetweenSlice * (userMatrixDimension3D - 1);
            getParam(FIELD_OF_VIEW_3D).setValue(fov3d);    // FOV ratio for display
            getParam(SPACING_BETWEEN_SLICE).setValue(spacingBetweenSlice);
        } else {
            sliceThickness = fov3d / userMatrixDimension3D;
            spacingBetweenSlice = 0;
            getParam(SPACING_BETWEEN_SLICE).setValue(spacingBetweenSlice);
            getParam(SLICE_THICKNESS).setValue(sliceThickness);
        }
    }

    @Override
    protected void getAcq4D() {
        // Avoid multi trigger time when  Multi echo or dynamic
        if (((ExtTrig) models.get("ExtTrig")).nb_trigger != 1 && (isDynamic)) {
            double tmp = ((ExtTrig) models.get("ExtTrig")).triggerTime.get(0);
            ((ExtTrig) models.get("ExtTrig")).triggerTime.clear();
            ((ExtTrig) models.get("ExtTrig")).triggerTime.add(tmp);
            ((ExtTrig) models.get("ExtTrig")).nb_trigger = 1;
        }
        if (isMultiplanar) {
            acqMatrixDimension4D = ((ExtTrig) models.get("ExtTrig")).nb_trigger * nb_dynamic_acquisition;
            userMatrixDimension4D = ((ExtTrig) models.get("ExtTrig")).nb_trigger * nb_dynamic_acquisition;
        } else {
            acqMatrixDimension4D = userMatrixDimension4D;
        }

        getParam(ACQUISITION_MATRIX_DIMENSION_4D).setValue(isKSCenterMode ? 1 : acqMatrixDimension4D);
        getParam(USER_MATRIX_DIMENSION_4D).setValue(userMatrixDimension4D);
    }

    @Override
    protected void getROGrad() throws Exception {
        gradReadout = Gradient.createGradient(getSequence(), Grad_amp_read_read, Time_rx, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
        if (isEnableRead && !gradReadout.calculateReadoutGradient(spectralWidth, getDouble(RESOLUTION_FREQUENCY) * acqMatrixDimension1D)) {
            double spectral_width_max = gradReadout.getSpectralWidth();
            if (getBoolean(SPECTRAL_WIDTH_OPT)) {
                notifyOutOfRangeParam(SPECTRAL_WIDTH, ((NumberParam) getParam(SPECTRAL_WIDTH)).getMinValue(), (spectral_width_max / (isFovDoubled ? 2 : 1)), "SPECTRAL_WIDTH too high for the readout gradient");
            } else {
                notifyOutOfRangeParam(SPECTRAL_WIDTH_PER_PIXEL, ((NumberParam) getParam(SPECTRAL_WIDTH_PER_PIXEL)).getMinValue(), (spectral_width_max / acqMatrixDimension1D), "SPECTRAL_WIDTH too high for the readout gradient");
            }
            spectralWidth = spectral_width_max;
        }
        gradReadout.applyAmplitude(Order.LoopB);
        set(Spectral_width, spectralWidth);
    }

    @Override
    protected void getRx() {
        pulseRX = RFPulse.createRFPulse(getSequence(), Time_rx, Rx_freq_offset, Rx_phase);
        pulseRX.setFrequencyOffsetReadout(gradReadout, off_center_distance_1D);

        //fill the OFF_CENTER_FIELD_OF_VIEW_EFF User Parameter
        ArrayList<Number> off_center_distanceList = new ArrayList<>();
        off_center_distanceList.add(off_center_distance_1D);
        off_center_distanceList.add(0);
        off_center_distanceList.add(0);

        getParam(OFF_CENTER_FIELD_OF_VIEW_EFF).setValue(off_center_distanceList);
    }

    @Override
    protected void getPhaseCyc() {
        //--------------------------------------------------------------------------------------
        //  calculate RF_SPOILING
        //--------------------------------------------------------------------------------------
        pulseRX.setPhase(0.0);

        RFPulse pulseRFSpoiler = RFPulse.createRFPulse(getSequence(), SP.Time_rf_spoiling, SP.FreqOffset_RFSpoiling);
        pulseRFSpoiler.setFrequencyOffsetForPhaseShift(is_rf_spoiling ? 117.0 : 0.0);

        // ----------------------------------------------------------------------------------------------
        // modify RX PHASE TABLE to handle OFF CENTER FOV 2D in both cases or PHASE CYCLING
        // ----------------------------------------------------------------------------------------------
        set(Rx_phase, 0);
    }

    @Override
    protected void getPrephaseGrad() {
        // -----------------------------------------------
        // calculate gradient equivalent rise time
        // -----------------------------------------------
        grad_shape_rise_time = clcGradEqRiseTime(Grad_shape_rise_up, Grad_shape_rise_down, grad_rise_time);
        set(Time_grad_phase_top, GRADIENT_PHASE_APPLICATION_TIME);

        // pre-calculate READ_prephasing max area
        gradReadPrep = Gradient.createGradient(getSequence(), Grad_amp_read_prep, Time_grad_phase_top, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
        if (isEnableRead) {
            gradReadPrep.refocalizeGradient(gradReadout, getDouble(PREPHASING_READ_GRADIENT_RATIO));
        }

        // pre-calculate SLICE_refocusing
        gradSliceRefPhase3D = Gradient.createGradient(getSequence(), Grad_amp_phase_3D_prep, Time_grad_phase_top, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
        if (isEnableSlice) {
            gradSliceRefPhase3D.refocalizeGradient(gradSlice, getDouble(SLICE_REFOCUSING_GRADIENT_RATIO));
        }
    }

    @Override
    protected void getPEGrad() {
        // pre-calculate PHASE_2D
        gradPhase2D = Gradient.createGradient(getSequence(), Grad_amp_phase_2D_prep, Time_grad_phase_top, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
        if (isEnablePhase) {
            gradPhase2D.preparePhaseEncoding(acqMatrixDimension2D, fovPhase, is_k_s_centred);
            gradPhase2D.reoderPhaseEncoding(plugin, 1, getInt(ACQUISITION_MATRIX_DIMENSION_2D), acqMatrixDimension1D);
        }

        // pre-calculate PHASE_3D
        if (!isMultiplanar && isEnablePhase3D) {
            gradSliceRefPhase3D.preparePhaseEncoding(acqMatrixDimension3D, slice_thickness_excitation, is_k_s_centred);
            gradSliceRefPhase3D.reoderPhaseEncoding3D(plugin, acqMatrixDimension3D);
        }
    }

    @Override
    protected void getGradOpt() {
        double grad_phase_application_time = getDouble(GRADIENT_PHASE_APPLICATION_TIME);
        double grad_area_sequence_max = 100 * (grad_phase_application_time + grad_shape_rise_time);
        double grad_area_max = Math.max(gradReadPrep.getTotalAbsArea(), Math.max(gradSliceRefPhase3D.getTotalAbsArea(), gradPhase2D.getTotalAbsArea()));            // calculate the maximum gradient aera between SLICE REFOC & READ PREPHASING
        if (grad_area_max > grad_area_sequence_max) {
            double grad_phase_application_time_min = ceilToSubDecimal(grad_area_max / 100.0 - grad_shape_rise_time, 6);
            notifyOutOfRangeParam(GRADIENT_PHASE_APPLICATION_TIME, grad_phase_application_time_min, ((NumberParam) getParam(GRADIENT_PHASE_APPLICATION_TIME)).getMaxValue(), "Gradient application time too short to reach this pixel dimension");
            grad_phase_application_time = grad_phase_application_time_min;
            set(Time_grad_phase_top, grad_phase_application_time);
            gradPhase2D.rePrepare();
            gradSliceRefPhase3D.rePrepare();
            gradReadPrep.rePrepare();
        }

        gradSliceRefPhase3D.applyAmplitude(Order.Three);
        gradPhase2D.applyAmplitude(Order.Two);
        gradReadPrep.applyAmplitude();
    }

    @Override
    protected void getSpoilerGrad() {
        // -------------------------------------------------------------------------------------------------
        // calculate Phase 2D, 3D and Read REWINDING - SPOILER area, check Grad_Spoil < GMAX
        // -------------------------------------------------------------------------------------------------

        // timing : grad_phase_application_time must be < grad_spoiler_application_time if rewinding
        //  boolean is_grad_rewinding = getBoolean(GRADIENT_ENABLE_REWINDING);// get slice refocussing ratio
        double grad_phase_application_time = getDouble(GRADIENT_PHASE_APPLICATION_TIME);
        double grad_spoiler_application_time = getDouble(GRADIENT_SPOILER_TIME);
        if (is_grad_rewinding && grad_phase_application_time > grad_spoiler_application_time) {
            notifyOutOfRangeParam(GRADIENT_SPOILER_TIME, grad_phase_application_time, ((NumberParam) getParam(GRADIENT_SPOILER_TIME)).getMaxValue(), "Gradient Spoiler top time must be longer than Phase Application Time");
            grad_spoiler_application_time = grad_phase_application_time;
        }
        set(Time_grad_spoiler_top, grad_spoiler_application_time);

        Gradient gradSliceSpoiler = Gradient.createGradient(getSequence(), Grad_amp_spoiler_slice, Time_grad_spoiler_top, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
        Gradient gradPhaseSpoiler = Gradient.createGradient(getSequence(), Grad_amp_spoiler_phase, Time_grad_spoiler_top, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
        Gradient gradReadSpoiler = Gradient.createGradient(getSequence(), Grad_amp_spoiler_read, Time_grad_spoiler_top, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);

        // Rewinding :
        if (is_grad_rewinding) {
            if (isEnablePhase3D)
                gradSliceSpoiler.refocalizePhaseEncodingGradient(gradSliceRefPhase3D);
            if (isEnableSlice)
                gradSliceSpoiler.refocalizeGradient(gradSlice, 1 - getDouble(SLICE_REFOCUSING_GRADIENT_RATIO));
            if (isEnablePhase)
                gradPhaseSpoiler.refocalizePhaseEncodingGradient(gradPhase2D);
            if (isEnableRead)
                gradReadSpoiler.refocalizeReadoutGradient(gradReadout, 1 - getDouble(PREPHASING_READ_GRADIENT_RATIO));
        }
        // Spoiler :
        //    List<Double> grad_amp_spoiler_sl_ph_re = getListDouble(grad_amp_spoiler_sl_ph_re);
        if (getBoolean(GRADIENT_ENABLE_SPOILER)) {
            double spoilerAmp = grad_amp_spoiler_sl_ph_re.get(0);
            if (!gradSliceSpoiler.addSpoiler(spoilerAmp))
                grad_amp_spoiler_sl_ph_re.set(0, spoilerAmp - gradSliceSpoiler.getSpoilerExcess());

            spoilerAmp = grad_amp_spoiler_sl_ph_re.get(1);
            if (!gradPhaseSpoiler.addSpoiler(spoilerAmp))
                grad_amp_spoiler_sl_ph_re.set(1, spoilerAmp - gradPhaseSpoiler.getSpoilerExcess());

            spoilerAmp = grad_amp_spoiler_sl_ph_re.get(2);
            if (!gradReadSpoiler.addSpoiler(spoilerAmp))
                grad_amp_spoiler_sl_ph_re.set(2, spoilerAmp - gradReadSpoiler.getSpoilerExcess());
        }
        gradPhaseSpoiler.applyAmplitude();
        gradSliceSpoiler.applyAmplitude();
        gradReadSpoiler.applyAmplitude();
    }

    @Override
    protected void getUPDisp() {
        this.getParam(TX_ATT).setValue(pulseTX.getAttParamValue());            // display PULSE_ATT
        this.getParam(TX_AMP).setValue(pulseTX.getAmp());     // display 90° amplitude
        this.getParam(TX_AMP_90).setValue(pulseTX.getAmp90());     // display 90° amplitude
        this.getParam(TX_AMP_180).setValue(pulseTX.getAmp180());   // display 180° amplitude
    }

    @Override
    protected void getTimeandDelay() throws Exception {
    }

    @Override
    protected void getTR() {
    }

    @Override
    protected void getAcqTime() {
        //----------------------------------------------------------------------
        // DYNAMIC SEQUENCEprepSeqTiming
        // Calculate frame acquisition time
        // Calculate delay between 4D acquisition
        //----------------------------------------------------------------------
        double frame_acquisition_time = nb_scan_1d * nb_scan_3d * nb_scan_2d * tr;
        double time_between_frames_min = ceilToSubDecimal(frame_acquisition_time + minInstructionDelay + min_flush_delay, 1);
        time_between_frames = time_between_frames_min;
        double interval_between_frames_delay = min_flush_delay;

        if (isDynamic && (acqMatrixDimension4D > 1)) {
            //Dynamic Sequence
            time_between_frames = getDouble(DYN_TIME_BTW_FRAMES);
            if (isDynamicMinTime) {
                time_between_frames = time_between_frames_min;
                getParam(DYN_TIME_BTW_FRAMES).setValue(time_between_frames_min);
            } else if (time_between_frames < (time_between_frames_min)) {
                notifyOutOfRangeParam(DYN_TIME_BTW_FRAMES, time_between_frames_min, ((NumberParam) getParam(DYN_TIME_BTW_FRAMES)).getMaxValue(), "Minimum frame acquisition time ");
                time_between_frames = time_between_frames_min;
            }
            interval_between_frames_delay = Math.max(time_between_frames - frame_acquisition_time, min_flush_delay);
        }
        set(Time_btw_dyn_frames, interval_between_frames_delay);
        // ------------------------------------------------------------------
        // Total Acquisition Time
        // ------------------------------------------------------------------
        double total_acquisition_time;
        if (!isMultiplanar) {
            total_acquisition_time = time_between_frames * Math.ceil(acqMatrixDimension4D / nb_interleaved_slice) + tr * nb_preScan;
        } else {
            total_acquisition_time = time_between_frames * nb_dynamic_acquisition + tr * nb_preScan;
        }
        getParam(SEQUENCE_TIME).setValue(total_acquisition_time);

    }

    @Override
    protected void getMultiParaList() {
        double frame_acquisition_time = nb_scan_1d * nb_scan_3d * nb_scan_2d * tr;

        int number_of_MultiSeries = 1;
        double time_between_MultiSeries = 0;
        ArrayList<Number> multiseries_valuesList = new ArrayList<>();
        String multiseries_parametername = "";

        if (models.get("ExtTrig").isEnabled() && ((ExtTrig) models.get("ExtTrig")).nb_trigger != 1) {
            number_of_MultiSeries = ((ExtTrig) models.get("ExtTrig")).nb_trigger;
            time_between_MultiSeries = frame_acquisition_time;
            multiseries_parametername = "TRIGGER DELAY";
            for (int i = 0; i < number_of_MultiSeries; i++) {
                double multiseries_value = Math.round(((ExtTrig) models.get("ExtTrig")).triggerTime.get(i) * 1e5) / 1e2;
                multiseries_valuesList.add(multiseries_value);
            }
        }
        getParam(MULTISERIES_PARAMETER_VALUE).setValue(multiseries_valuesList);
        getParam(MULTISERIES_PARAMETER_NAME).setValue(multiseries_parametername);

        ArrayList<Number> acquisition_timesList = new ArrayList<>();
        double acqusition_time;
        for (int i = 0; i < nb_dynamic_acquisition; i++) {
            for (int j = 0; j < number_of_MultiSeries; j++) {
                acqusition_time = roundToDecimal((i * time_between_frames + j * time_between_MultiSeries), 3);
                acquisition_timesList.add(acqusition_time);
            }
        }
        getParam(ACQUISITION_TIME_OFFSET).setValue(acquisition_timesList);
    }

    protected void getGradRiseTime() {
        double min_rise_time_factor = getDouble(MIN_RISE_TIME_FACTOR);

        double min_rise_time_sinus = GradientMath.getShortestRiseTime(100.0) * Math.PI / 2 * 100 / min_rise_time_factor;
        if (grad_rise_time < min_rise_time_sinus) {
            double new_grad_rise_time = ceilToSubDecimal(min_rise_time_sinus, 5);
            notifyOutOfRangeParam(GRADIENT_RISE_TIME, new_grad_rise_time, ((NumberParam) getParam(GRADIENT_RISE_TIME)).getMaxValue(), "Gradient ramp time too short ");
            grad_rise_time = new_grad_rise_time;
        }
        set(Time_grad_ramp, grad_rise_time);
    }

}
