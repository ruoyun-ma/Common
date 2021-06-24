package kernel;

import rs2d.spinlab.instrument.util.GradientMath;
import rs2d.spinlab.sequence.SequenceTool;
import rs2d.spinlab.sequence.table.Table;
import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.*;
import rs2d.spinlab.tools.table.Order;

import java.util.*;

import static java.util.Arrays.asList;

import common.*;
import model.*;

import static common.CommonSP.*;
import static common.CommonUP.*;
import static kernel.KernelSE.SP.*;
import static kernel.KernelSE.UP.*;

/**
 * Abstract Class KernelGE
 * prep Kernel functions of Gradient Echo Sequences
 * V1.0- 2021-4-26 XG
 * <p>
 * Keyhole/ TOF/ Flow/ Elliptical Sampling etc. are kept in each individual sequence
 * <p>
 */

public abstract class KernelSE extends SeqPrep {
    protected String sequenceVersion = "Version KernelSE 1.0";
    protected boolean CameleonVersion105 = false;
    protected boolean is_FSE_vs_MultiEcho;

    protected double min_flush_delay;
    protected double slice_thickness_excitation;
    protected int echoEffective;
    protected double txLength90;
    protected double txLength180;
//    protected List<Double> grad_amp_spoiler_sl_ph_re;

    //Dynamic
    protected boolean isDynamic;
    protected int nb_dynamic_acquisition;
    protected boolean isDynamicMinTime;
    protected double time_between_frames;

    protected boolean is_grad_spoiler;
    protected boolean is_k_s_centred;

    protected double grad_shape_rise_time;
    protected RFPulse pulseRX;
    protected Gradient gradPhase2D;
    protected Gradient gradReadPrep;
    protected Gradient gradReadout;
    protected Gradient gradSliceRefPhase3D;

    protected boolean isFSETrain1D;
    protected String transformplugin;
    protected double[] TE_TR_lim = {0, 100000, 0, 100000}; // limite {TEmin TEmax, TRmin TRmax}


    protected enum UP implements GeneratorParamEnum {
        TOOL_FSE_TRAIN_1D,
        KS_CENTERED,
        DYNAMIC_SEQUENCE,
        DYN_NUMBER_OF_ACQUISITION,
        DYNAMIC_MIN_TIME,
        DYN_TIME_BTW_FRAMES,
        GRADIENT_ENABLE_SPOILER,
        GRADIENT_ENABLE_SLICE_CRUSH,
        GRADIENT_ENABLE_SLICE_CRUSHER,
        GRAD_AMP_SPOILER_SL_PH_RE,
        SLICE_REFOCUSING_GRADIENT_RATIO,
        NUMBER_OF_INTERLEAVED_SLICE,
        GRADIENT_PHASE_APPLICATION_TIME,
        GRADIENT_SPOILER_TIME,
        GRADIENT_CRUSHER_TOP_TIME,
        TRIGGER_TIME,
        IMAGE_CONTRAST,
        LIM_T1_WEIGHTED,
        LIM_T2_WEIGHTED,
        LIM_PD_WEIGHTED,
        PHASE_CYCLING,
        SLICE_THICKNESS_180_WIDER,
        GRADIENT_READ_OFFSET,
        GRADIENT_READ_PREPHASING_APPLICATION_TIME,
        GRADIENT_AREA_CRUSHER,
        GRADIENT_AREA_CRUSHER_PI,
        GRADIENT_AMP_CRUSHER,
        GRADIENT_CRUSHER_END_TOP_TIME,
        GRADIENT_AMP_SPOILER,
        GRADIENT_CRUSHER_READ_TOP_TIME,
        ECHO_EFFECTIVE,
        SE_TYPE,

        ;

        @Override
        public Param build() {
            return null;
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum {
        Grad_enable_read,
        Grad_enable_phase_2D,
        Grad_enable_phase_3D,
        Grad_enable_slice,
        Grad_enable_slice_180,
        Grad_enable_slice_crush,
        Grad_enable_spoiler,
        Grad_amp_read_read,
        Grad_amp_slice_180,
        FreqOffset_tx_prep,
        FreqOffset_tx_comp,
        FreqOffset_tx_prep_180,
        Frequency_offset_init,
        Time_btw_dyn_frames,

        Grad_amp_spoiler_slice,
        Grad_amp_spoiler_phase,
        Grad_amp_spoiler_read,
        Grad_amp_slice_refoc,
        Grad_amp_slice_crush,
        Time_grad_spoiler_top,
        Time_grad_crusher_top,
        Time_grad_read_prep_top,
        Time_grad_read_crusher,

        //Below UP maybe removed in next version
        ;
    }

    public KernelSE() {
        super();
        addUserParams();
    }

    @Override
    public void init() {
        super.init();
        ((TextParam) getParam(TX_SHAPE_90)).setSuggestedValues(tx_shape);
        ((TextParam) getParam(TX_SHAPE_90)).setRestrictedToSuggested(true);
        ((TextParam) getParam(TX_SHAPE_180)).setSuggestedValues(tx_shape);
        ((TextParam) getParam(TX_SHAPE_180)).setRestrictedToSuggested(true);
        //TRANSFORM PLUGIN
        TextParam transformPlugin = getParam(CommonUP.TRANSFORM_PLUGIN);
        transformPlugin.setSuggestedValues(asList("Centered2DRot",
                "Bordered2D",
                "Sequential4D",
                "Sequential2D",
                "FSE_TRAIN_1D",
                "Sequential2DInterleaved"));
        transformPlugin.setRestrictedToSuggested(true);
    }

    // ==============================
// -----   GENERATE
// ==============================
    @Override
    public void initUserParam() {
        super.initUserParam();
        getParam(SEQUENCE_VERSION).setValue(sequenceVersion);

        txLength90 = getDouble(TX_LENGTH_90);
        txLength180 = getDouble(TX_LENGTH_180);

        is_FSE_vs_MultiEcho = ("FSE".equalsIgnoreCase((String) (getParam(SE_TYPE).getValue())));
        nb_echo_4D = !is_FSE_vs_MultiEcho ? echoTrainLength : 1; // didn't support vfl mode

        isDynamic = getBoolean(DYNAMIC_SEQUENCE);
        nb_dynamic_acquisition = isDynamic ? getInt(DYN_NUMBER_OF_ACQUISITION) : 1;
        isDynamic = isDynamic && (nb_dynamic_acquisition > 1);
        isDynamicMinTime = getBoolean(DYNAMIC_MIN_TIME);

        is_grad_spoiler = getBoolean(GRADIENT_ENABLE_SPOILER);// get slice refocussing ratio
        is_k_s_centred = getBoolean(KS_CENTERED);

        isFSETrain1D = getBoolean(TOOL_FSE_TRAIN_1D);
        isFSETrain1D = !isKSCenterMode && isFSETrain1D;
        getParam(TOOL_FSE_TRAIN_1D).setValue(isFSETrain1D);
        isEnablePhase3D = !isKSCenterMode && !isFSETrain1D && isEnablePhase3D;
        isEnablePhase = !isKSCenterMode && !isFSETrain1D && isEnablePhase;
    }

    //--------------------------------------------------------------------------------------
    // ini functions for beforeRouting() and initAfterRouting()
    //--------------------------------------------------------------------------------------
    @Override
    protected void iniModels() {
        setModels(new ArrayList<>(Arrays.asList(new ExtTrig())), this);
    }

    @Override
    protected void iniTransformPlugin() throws Exception {
        echoEffective = Math.min(echoTrainLength, echoEffective);

        //  get limits for the image contrast
        switch (getText(IMAGE_CONTRAST)) {
            case "T1-weighted": // Short TE, Short TR
                ListNumberParam T1_contrast_TE_TR_lim = getParam(LIM_T1_WEIGHTED);
                TE_TR_lim[1] = T1_contrast_TE_TR_lim.getValue().get(0).doubleValue(); //TEmax
                TE_TR_lim[2] = T1_contrast_TE_TR_lim.getValue().get(1).doubleValue(); //TRmin
                TE_TR_lim[3] = T1_contrast_TE_TR_lim.getValue().get(2).doubleValue(); //TRmax
                break;
            case "PD-weighted": // Short TE, Long TR
                ListNumberParam PD_contrast_TE_TR_lim = getParam(LIM_PD_WEIGHTED);
                TE_TR_lim[1] = PD_contrast_TE_TR_lim.getValue().get(0).doubleValue(); //TEmax
                TE_TR_lim[2] = PD_contrast_TE_TR_lim.getValue().get(1).doubleValue(); //TRmin
                break;
            case "T2-weighted": // Long TE, Long TR
                ListNumberParam T2_contrast_TE_TR_lim = getParam(LIM_T2_WEIGHTED);
                TE_TR_lim[0] = T2_contrast_TE_TR_lim.getValue().get(0).doubleValue(); //TEmin
                TE_TR_lim[2] = T2_contrast_TE_TR_lim.getValue().get(1).doubleValue(); //TRmin
                break;
            default: // Custom
                break;
        }

        ///// transform plugin
        if (isKSCenterMode) {
            transformplugin = "Sequential2D";
        } else if (isFSETrain1D) {
            transformplugin = "FSE_TRAIN_1D";
        } else if (is_FSE_vs_MultiEcho) {
            switch (getText(IMAGE_CONTRAST)) {
                case "T1-weighted":
                case "PD-weighted":
                    transformplugin = "Centered2DRot"; // first echo should be the significant echo
                    echoEffective = 1;
                    break;
                case "T2-weighted":
                    transformplugin = "Centered2DRot"; // first echo should be the significant echo
                    echoEffective = Math.max((int) Math.ceil(echoTrainLength / 2), echoEffective);
                    break;
                default: // Custom
                    switch (transformplugin) {
                        case "Centered2D": // not used anymore, replace by Centered2DRot
                        case "Sequential4D":
                            transformplugin = "Centered2DRot";
                            echoEffective = 1;
                            break;
                        case "Bordered2D": // not used anymore, replace by Centered2DRot
                            transformplugin = "Centered2DRot";
                            echoEffective = echoTrainLength;
                            break;
                        case "Sequential2D":
                            if (echoTrainLength != 1)// not isKSCenterMode anymore, restor FSE
                                transformplugin = "Centered2DRot";
                            break;
                        case "FSE_TRAIN_1D":// not isFSETrainTest anymore, restor FSE
                            transformplugin = "Centered2DRot";
                            break;
                        default:
                            break;
                    }
                    break;
            }
            te = echoEffective * echo_spacing;
        } else { // is_FSE_vs_MultiEcho = false  :  multi echo for T2 measurement
            transformplugin = "Sequential4D";
            getParam(IMAGE_CONTRAST).setValue("T2-weighted");
            echo_spacing = te;
        }

        getParam(ECHO_SPACING).setValue(echo_spacing);
        getParam(ECHO_EFFECTIVE).setValue(echoEffective);
        getParam(TRANSFORM_PLUGIN).setValue(transformplugin);

        plugin = getTransformPlugin();
        plugin.setParameters(new ArrayList<>(getUserParams()));
        plugin.getScanOrder(); //output traj. graphs for Elliptical3D plugin
    }

    @Override
    protected void iniScanLoop() {
        nb_scan_1d = nb_averages;
        nb_scan_2d = acqMatrixDimension2D;
        String updateDimension = SequenceTool.UpdateDimensionEnum.Disable.getText();
        if (!isMultiplanar) {
            nb_scan_3d = acqMatrixDimension3D;
        } else {
            nb_scan_3d = nb_shoot_3d;
        }

        //nb_scan_4d = models.get(ExtTrig.class).nb_trigger * nb_dynamic_acquisition;
        nb_scan_4d = models.get(ExtTrig.class).nb_trigger
                * models.get(InvRec.class).nb_inversionRecovery
                * nb_dynamic_acquisition;


        nb_planar_excitation = (isMultiplanar ? acqMatrixDimension3D : 1);

        if (isKSCenterMode) { // Do only the center of the k-space for auto RG
            nb_scan_1d = 1;
            nb_scan_2d = 2;
            nb_scan_3d = !isMultiplanar ? 1 : nb_scan_3d;
            nb_scan_4d = 1;
            getParam(ACQUISITION_MATRIX_DIMENSION_3D).setValue(!isMultiplanar ? 1 : acqMatrixDimension3D);
        }
        if (isFSETrain1D) { // Do only the center of the k-space for auto RG
            nb_scan_1d = 1;
            nb_scan_3d = !isMultiplanar ? 1 : nb_scan_3d;
            nb_scan_4d = 1;
            getParam(ACQUISITION_MATRIX_DIMENSION_1D).setValue(acqMatrixDimension1D * echoTrainLength);
            getParam(ACQUISITION_MATRIX_DIMENSION_2D).setValue(nb_scan_2d);
            getParam(ACQUISITION_MATRIX_DIMENSION_3D).setValue(!isMultiplanar ? 1 : acqMatrixDimension3D);
        }

        // set Nb_scan  Values
        set(Pre_scan, nb_preScan); // Do the prescan
        set(Nb_point, acqMatrixDimension1D);
        set(Nb_1d, nb_scan_1d);
        set(Nb_2d, nb_scan_2d);
        set(Nb_3d, nb_scan_3d);
        //set(Update_dimension, updateDimension);
        set(Nb_4d, nb_scan_4d);
        set(Nb_echo, echoTrainLength - 1);
        set(Nb_interleaved_slice, nb_interleaved_slice - 1);
    }

    @Override
    protected void iniSeqDisp() {
        String seqDescription = "SE_";
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
        if (echoTrainLength != 1) {
            seqDescription += "_ETL=" + echoTrainLength;
        }
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
        set(Grad_enable_phase_3D, (!isMultiplanar && isEnablePhase3D));
        set(Grad_enable_slice, isEnableSlice);
        set(Grad_enable_slice_180, isEnableSlice);

        if (hasParam(GRADIENT_ENABLE_SLICE_CRUSHER))
            set(Grad_enable_slice_crush, GRADIENT_ENABLE_SLICE_CRUSHER);
        else {
            set(Grad_enable_slice_crush, GRADIENT_ENABLE_SLICE_CRUSH);
        }
        set(Grad_enable_spoiler, GRADIENT_ENABLE_SPOILER);
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
        set(Time_tx, TX_LENGTH_90);     // set RF pulse length to sequence
        set(Time_tx_180, TX_LENGTH_180);   // set 180° RF pulse length to sequence
        pulseTX = RFPulse.createRFPulse(getSequence(), Tx_att, Tx_amp, Tx_phase, Time_tx, Tx_shape, Tx_shape_phase, Tx_freq_offset);
        pulseTX180 = RFPulse.createRFPulse(getSequence(), Tx_att, Tx_amp_180, Tx_phase_180, Time_tx_180, Tx_shape_180, Tx_shape_phase_180, Tx_freq_offset_180);

        if (getSequence().getPublicTable(Tx_att_offset.name()) != null) {
            pulseTX.createAttOffset(getSequence(), Tx_att_offset);
        }
        if (getSequence().getPublicTable(Tx_att_offset_180.name()) != null) {
            pulseTX180.createAttOffset(getSequence(), Tx_att_offset_180);
        }
        double flip_angle = getDouble(FLIP_ANGLE);

        // -----------------------------------------------
        // Calculation RF pulse parameters  2/4 : Shape
        // -----------------------------------------------
        pulseTX.setShape((getText(TX_SHAPE_90)), nb_shape_points, "90 degree");
        pulseTX180.setShape((getText(TX_SHAPE_180)), nb_shape_points, "Refocusing (spin-echo)");

        // -----------------------------------------------
        // Calculation RF pulse parameters  3/4 : RF pulse & attenuation
        // -----------------------------------------------
        if (getBoolean(TX_AMP_ATT_AUTO)) {
            if (!pulseTX.checkPower(flip_angle, observeFrequency, nucleus)) {
                notifyOutOfRangeParam(TX_LENGTH_90, pulseTX.getPulseDuration(), ((NumberParam) getParam(TX_LENGTH_90)).getMaxValue(), "Pulse length too short to reach RF power with this pulse shape");
                txLength90 = pulseTX.getPulseDuration();
            }
            if (!pulseTX180.checkPower(180, observeFrequency, nucleus)) {
                notifyOutOfRangeParam(TX_LENGTH_180, pulseTX180.getPulseDuration(), ((NumberParam) getParam(TX_LENGTH_180)).getMaxValue(), "Pulse length too short to reach RF power with this pulse shape");
                txLength180 = pulseTX180.getPulseDuration();
            }
            pulseTX180.prepAtt(80, getListInt(TX_ROUTE));
            pulseTX.prepTxAmp(getListInt(TX_ROUTE));
            pulseTX180.prepTxAmp(getListInt(TX_ROUTE));
            rfPulses.add(pulseTX);
            rfPulses.add(pulseTX180);
            getUPDisp();
        } else {
            pulseTX.setAtt(getInt(TX_ATT));
            pulseTX.setAmp(getDouble(TX_AMP_90));
            pulseTX180.setAmp(getParam(TX_AMP_180));
        }

        // -----------------------------------------------
        // Calculation RF pulse parameters  4/4: bandwidth
        // -----------------------------------------------
        double tx_bandwidth_factor_90 = getTx_bandwidth_factor(TX_SHAPE_90, TX_BANDWIDTH_FACTOR, TX_BANDWIDTH_FACTOR_3D);
        double tx_bandwidth_90 = tx_bandwidth_factor_90 / txLength90;

        double tx_bandwidth_factor_180 = getTx_bandwidth_factor(TX_SHAPE_180, TX_BANDWIDTH_FACTOR, TX_BANDWIDTH_FACTOR_3D);
        double tx_bandwidth_180 = tx_bandwidth_factor_180 / txLength180;

        // ---------------------------------------------------------------------
        // calculate SLICE gradient amplitudes for RF pulses
        // ---------------------------------------------------------------------
        slice_thickness_excitation = (isMultiplanar ? sliceThickness : (sliceThickness * userMatrixDimension3D));
        boolean is_wide_180_slice_thickness = getBoolean(SLICE_THICKNESS_180_WIDER);
        double slice_thickness_excitation_180;
        if (((spacingBetweenSlice > slice_thickness_excitation) || !isMultiplanar) && is_wide_180_slice_thickness) {
            slice_thickness_excitation_180 = slice_thickness_excitation * (1 + 0.5);
        } else if ((spacingBetweenSlice != 0) && is_wide_180_slice_thickness) {
            slice_thickness_excitation_180 = slice_thickness_excitation + spacingBetweenSlice / 2.0;
        } else {
            slice_thickness_excitation_180 = slice_thickness_excitation;
        }

        gradSlice = Gradient.createGradient(getSequence(), Grad_amp_slice, Time_tx, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
        gradSlice180 = Gradient.createGradient(getSequence(), Grad_amp_slice_180, Time_tx_180, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);

        if (isEnableSlice) {
            boolean test90 = gradSlice.prepareSliceSelection(tx_bandwidth_90, slice_thickness_excitation);
            boolean test180 = gradSlice180.prepareSliceSelection(tx_bandwidth_180, slice_thickness_excitation_180);
            if (!test90 || !test180) {
                slice_thickness_excitation = Math.max(gradSlice.getSliceThickness(), gradSlice180.getSliceThickness());
                double slice_thickness_min = (isMultiplanar ? slice_thickness_excitation : (slice_thickness_excitation / userMatrixDimension3D));
                notifyOutOfRangeParam(SLICE_THICKNESS, slice_thickness_min, ((NumberParam) getParam(SLICE_THICKNESS)).getMaxValue(), "Pulse length too short to reach this slice thickness");
                sliceThickness = slice_thickness_min;
            }
        }
        gradSlice.applyAmplitude();
        gradSlice180.applyAmplitude(slice_thickness_excitation_180);

        // ------------------------------------------------------------------
        //calculate TX FREQUENCY offsets tables for slice positionning
        // ------------------------------------------------------------------
        if (isMultiplanar && nb_planar_excitation > 1 && isEnableSlice) {
            //MULTI-PLANAR case : calculation of frequency offset table
            pulseTX.prepareOffsetFreqMultiSlice(gradSlice, nb_planar_excitation, spacingBetweenSlice, off_center_distance_3D);
            pulseTX.reoderOffsetFreq(plugin, acqMatrixDimension1D * echoTrainLength, nb_interleaved_slice);
            pulseTX.setFrequencyOffset(nb_interleaved_slice != 1 ? Order.ThreeLoop : Order.Three);
            pulseTX180.prepareOffsetFreqMultiSlice(gradSlice180, nb_planar_excitation, spacingBetweenSlice, off_center_distance_3D);
            pulseTX180.reoderOffsetFreq(plugin, acqMatrixDimension1D * echoTrainLength, nb_interleaved_slice);
            pulseTX180.setFrequencyOffset(nb_interleaved_slice != 1 ? Order.ThreeLoop : Order.Three);
        } else {
            //3D CASE :
            pulseTX.prepareOffsetFreqMultiSlice(gradSlice, 1, 0, off_center_distance_3D);
            pulseTX.setFrequencyOffset(Order.Three);
            pulseTX180.prepareOffsetFreqMultiSlice(gradSlice180, 1, 0, off_center_distance_3D);
            pulseTX180.setFrequencyOffset(Order.Three);
        }

        // ------------------------------------------------------------------
        // calculate TX FREQUENCY offsets compensation
        // ------------------------------------------------------------------
        double grad_ratio_slice_refoc = isEnableSlice ? getDouble(SLICE_REFOCUSING_GRADIENT_RATIO) : 0.0;   // get slice refocussing ratio

        RFPulse pulseTXPrep = RFPulse.createRFPulse(getSequence(), Time_grad_ramp, FreqOffset_tx_prep);
        pulseTXPrep.setCompensationFrequencyOffset(pulseTX, grad_ratio_slice_refoc);

        RFPulse pulseTX180Prep = RFPulse.createRFPulse(getSequence(), Time_grad_ramp, FreqOffset_tx_prep_180);
        pulseTX180Prep.setCompensationFrequencyOffset(pulseTX180, grad_ratio_slice_refoc);

    }

    @Override
    protected void prepDicom() {
        // Set  ECHO_TIME
        if (!is_FSE_vs_MultiEcho) {
            ArrayList<Number> arrayListEcho = new ArrayList<>();
            for (int i = 0; i < acqMatrixDimension4D; i++) {
                arrayListEcho.add(te * i);
            }
            ListNumberParam list = new ListNumberParam(getParam(MriDefaultParams.ECHO_TIME), arrayListEcho);       // associate TE to images for DICOM export
            putVariableParameter(list, (4));
        }
    }

    //--------------------------------------------------------------------------------------
    // get functions
    //--------------------------------------------------------------------------------------
    @Override
    protected void getAcq2D() {
        super.getAcq2D();

        double partial_phase = getDouble(USER_PARTIAL_PHASE);
        double zero_filling_2D = (100 - partial_phase) / 100f;
        if (is_FSE_vs_MultiEcho) {
            int shoot = (int) Math.round((1 - zero_filling_2D) * userMatrixDimension2D / (2.0 * echoTrainLength));
            zero_filling_2D = 1.0 - shoot * (2.0 * echoTrainLength) / ((double) userMatrixDimension2D);
            if (zero_filling_2D < 0)
                zero_filling_2D = 1.0 - (shoot - 1) * (2.0 * echoTrainLength) / ((double) userMatrixDimension2D);
            partial_phase = (100 - zero_filling_2D * 100f);
        }
        getParam(USER_ZERO_FILLING_2D).setValue((100 - partial_phase));
        // -----------------------------------------------
        // 2nd D managment  ETL  //TODO:XG:Optimize the logic flow
        // -----------------------------------------------
        if (is_FSE_vs_MultiEcho) {
            echoTrainLength = getInferiorDivisorToGetModulusZero(echoTrainLength, acqMatrixDimension2D / 2);
            getParam(ECHO_TRAIN_LENGTH).setValue(echoTrainLength);
            nb_scan_2d = Math.floorDiv(acqMatrixDimension2D, echoTrainLength);
            if (transformplugin != null && !(transformplugin.equalsIgnoreCase("Sequential2D"))) {
                nb_scan_2d = floorEven(nb_scan_2d); // Centred or bordered skim need to be multiple of 2
            }
            acqMatrixDimension2D = nb_scan_2d * echoTrainLength;
        }
        userMatrixDimension2D = Math.max(userMatrixDimension2D, acqMatrixDimension2D);
        getParam(USER_MATRIX_DIMENSION_2D).setValue(userMatrixDimension2D);
    }

    @Override
    protected void getAcq3D() {
        // 3D ZERO FILLING
        double partial_slice;
        if (isMultiplanar) {
            getParam(USER_PARTIAL_SLICE).setValue(100);
            partial_slice = 100;
        } else {
            partial_slice = getDouble(USER_PARTIAL_SLICE);
        }
        getParam(USER_ZERO_FILLING_3D).setValue((100 - partial_slice));

        //Calculate the number of k-space lines acquired in the 3rd Dimension : acqMatrixDimension3D
        if (!isMultiplanar) {
            acqMatrixDimension3D = floorEven(partial_slice / 100f * userMatrixDimension3D);
            acqMatrixDimension3D = (acqMatrixDimension3D < 4) && isEnablePhase3D ? 4 : acqMatrixDimension3D;
            userMatrixDimension3D = Math.max(userMatrixDimension3D, acqMatrixDimension3D);
        } else {
            acqMatrixDimension3D = userMatrixDimension3D;
        }

        if (!isMultiplanar) {
            nb_shoot_3d = 0;
            nb_interleaved_slice = 1;
            nb_planar_excitation = 1;
        } else {
            nb_shoot_3d = getInferiorDivisorToGetModulusZero(nb_shoot_3d, acqMatrixDimension3D);
            nb_interleaved_slice = (int) Math.ceil((acqMatrixDimension3D / (double) nb_shoot_3d));
            nb_planar_excitation = userMatrixDimension3D;
        }

        getParam(NUMBER_OF_SHOOT_3D).setValue(nb_shoot_3d);
        //getParam(NUMBER_OF_INTERLEAVED_SLICE).setValue(isMultiplanar ? nb_interleaved_slice : 0);
        getParam(NUMBER_OF_INTERLEAVED_SLICE).setValue(nb_interleaved_slice); //XG
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
        if (models.get(ExtTrig.class).nb_trigger != 1 &&
                (models.get(InvRec.class).nb_inversionRecovery != 1 || isDynamic || (nb_echo_4D != 1))) {
            double tmp = models.get(ExtTrig.class).triggerTime.get(0);
            models.get(ExtTrig.class).triggerTime.clear();
            models.get(ExtTrig.class).triggerTime.add(tmp);
            models.get(ExtTrig.class).nb_trigger = 1;
        }

        userMatrixDimension4D = models.get(ExtTrig.class).nb_trigger
                * models.get(InvRec.class).nb_inversionRecovery
                * nb_dynamic_acquisition
                * nb_echo_4D;

        acqMatrixDimension4D = userMatrixDimension4D;

        getParam(ACQUISITION_MATRIX_DIMENSION_4D).setValue(isKSCenterMode ? 1 : acqMatrixDimension4D);
        getParam(USER_MATRIX_DIMENSION_4D).setValue(userMatrixDimension4D);
    }

    @Override
    protected void getROGrad() throws Exception {
        double grad_crusher_read_time = getDouble(GRADIENT_CRUSHER_READ_TOP_TIME);
        set(Time_grad_read_crusher, grad_crusher_read_time);

        gradReadout = Gradient5Event.createGradient(getSequence(), Grad_amp_read_read, Time_grad_read_crusher, Time_rx, Time_grad_read_crusher, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
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
        // modify RX PHASE TABLE to handle OFF CENTER FOV 2D in both cases or PHASE CYCLING
        //--------------------------------------------------------------------------------------
        pulseRX.setPhase(0.0);

        boolean is_phase_cycling = getBoolean(PHASE_CYCLING);
        if (is_phase_cycling) { // to do modify the phase cycling
            pulseTX.setPhase(nb_averages != 1 ? Order.One : Order.Two, 0.0, 180.0);
            pulseRX.setPhase(nb_averages != 1 ? Order.One : Order.Two, 0.0, 180.0);
        } else {
            pulseTX.setPhase(0.0);
            pulseRX.setPhase(0.0);
        }
        pulseTX180.setPhase(90);
    }

    @Override
    protected void getPrephaseGrad() {
        // -----------------------------------------------
        // calculate gradient equivalent rise time
        // -----------------------------------------------
        grad_shape_rise_time = clcGradEqRiseTime(Grad_shape_rise_up, Grad_shape_rise_down, grad_rise_time);

        double grad_read_prep_application_time;
        double grad_read_prep_offset = getDouble(GRADIENT_READ_OFFSET);
        //set(Time_grad_read_prep_top, grad_read_prep_application_time);

        // calc grad_read_prep_application_time by hand  TODO: Merge into Gradient.class
        double tmp_StaticArea = getDouble(PREPHASING_READ_GRADIENT_RATIO) * gradReadout.getStaticArea();
        grad_read_prep_application_time = tmp_StaticArea / gradReadout.getAmplitude() - grad_shape_rise_time;

        getParam(GRADIENT_READ_PREPHASING_APPLICATION_TIME).setValue(grad_read_prep_application_time);
        set(Time_grad_read_prep_top, grad_read_prep_application_time);
        set(Grad_amp_read_prep, gradReadout.getAmplitude());

        // pre-calculate READ_prephasing max area
        Gradient gradReadPrep = Gradient.createGradient(getSequence(), Grad_amp_read_prep, Time_grad_read_prep_top, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
        if (isEnableRead)
            gradReadPrep.refocalizeGradient(gradReadout, -getDouble(PREPHASING_READ_GRADIENT_RATIO));
        if (!gradReadPrep.addSpoiler(grad_read_prep_offset))
            getParam(GRADIENT_READ_OFFSET).setValue(grad_read_prep_offset - gradReadPrep.getSpoilerExcess());   // display observation time

        // pre-calculate SLICE_refocusing
        double grad_ratio_slice_refoc = 0.5;   // get slice refocussing ratio
        this.getParam(SLICE_REFOCUSING_GRADIENT_RATIO).setValue(grad_ratio_slice_refoc);   // display 180° amplitude
        Gradient gradSliceRef = Gradient.createGradient(getSequence(), Grad_amp_slice_refoc, Time_grad_read_prep_top, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
        if (isEnableSlice) {
            gradSliceRef.refocalizeGradient(gradSlice, grad_ratio_slice_refoc);
        }

        // Check if enougth time for 2D_PHASE, 3D_PHASE SLICE_REF or READ_PREP
        double grad_area_sequence_max = 100 * (grad_read_prep_application_time + grad_shape_rise_time);
        double grad_area_max = Math.max(gradReadPrep.getTotalArea(), gradSliceRef.getTotalArea());            // calculate the maximum gradient aera between SLICE REFOC & READ PREPHASING
        if (grad_area_max > grad_area_sequence_max) {
            double grad_read_prep_application_time_min = ceilToSubDecimal(grad_area_max / 100.0 - grad_shape_rise_time, 5);
            notifyOutOfRangeParam(GRADIENT_READ_PREPHASING_APPLICATION_TIME, grad_read_prep_application_time_min, ((NumberParam) getParam(GRADIENT_READ_PREPHASING_APPLICATION_TIME)).getMaxValue(), "Gradient application time too short to reach this pixel dimension");
            grad_read_prep_application_time = grad_read_prep_application_time_min;
            set(Time_grad_read_prep_top, grad_read_prep_application_time);
            gradSliceRef.rePrepare();
            gradReadPrep.rePrepare();
        }
        gradSliceRef.applyAmplitude();
        gradReadPrep.applyAmplitude();
    }

    @Override
    protected void getPEGrad() {
        // -------------------------------------------------------------------------------------------------
        // calculate PHASE_3D  & PHASE_2D
        // -------------------------------------------------------------------------------------------------
        double grad_phase_application_time = getDouble(GRADIENT_PHASE_APPLICATION_TIME);
        set(Time_grad_phase_top, grad_phase_application_time);

        // pre-calculate PHASE_3D
        gradSliceRefPhase3D = Gradient.createGradient(getSequence(), Grad_amp_phase_3D_prep, Time_grad_phase_top, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
        if (!isMultiplanar && isEnablePhase3D)
            gradSliceRefPhase3D.preparePhaseEncodingForCheck(acqMatrixDimension3D, acqMatrixDimension3D, slice_thickness_excitation, getBoolean(KS_CENTERED));

        // pre-calculate PHASE_2D
        gradPhase2D = Gradient.createGradient(getSequence(), Grad_amp_phase_2D_prep, Time_grad_phase_top, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
        gradPhase2D.preparePhaseEncodingForCheck(acqMatrixDimension2D, acqMatrixDimension2D, fovPhase, getBoolean(KS_CENTERED));
        if (is_FSE_vs_MultiEcho && echoTrainLength != 1 && !isKSCenterMode && !isFSETrain1D) {
            gradPhase2D.reoderPhaseEncoding(plugin, echoTrainLength, acqMatrixDimension2D, acqMatrixDimension1D);
        }

        // Check if enougth time for 2D_PHASE, 3D_PHASE SLICE_REF or READ_PREP
        double grad_area_sequence_max = 100 * (grad_phase_application_time + grad_shape_rise_time);
        double grad_area_max = Math.max(gradSliceRefPhase3D.getTotalArea(), gradPhase2D.getTotalArea());            // calculate the maximum gradient aera between SLICE REFOC & READ PREPHASING

        if (grad_area_max > grad_area_sequence_max) {
            double grad_phase_application_time_min = ceilToSubDecimal(grad_area_max / 100.0 - grad_shape_rise_time, 5);
            notifyOutOfRangeParam(GRADIENT_PHASE_APPLICATION_TIME, grad_phase_application_time_min, ((NumberParam) getParam(GRADIENT_PHASE_APPLICATION_TIME)).getMaxValue(), "Gradient application time too short to reach this pixel dimension");
            grad_phase_application_time = grad_phase_application_time_min;
            set(Time_grad_phase_top, grad_phase_application_time);
            gradPhase2D.rePrepare();
            gradSliceRefPhase3D.rePrepare();
        }
        gradSliceRefPhase3D.applyAmplitude(Order.Three);
        gradPhase2D.applyAmplitude((is_FSE_vs_MultiEcho && echoTrainLength != 1) ? Order.TwoLoopB : Order.Two);

        // Reorder Phase Encoding
        Gradient gradSlicePhase3D_comp = Gradient.createGradient(getSequence(), Grad_amp_phase_3D_comp, Time_grad_phase_top, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
        Gradient gradPhase2D_comp = Gradient.createGradient(getSequence(), Grad_amp_phase_2D_comp, Time_grad_phase_top, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);

        if (!isMultiplanar && isEnablePhase3D)
            gradSlicePhase3D_comp.refocalizePhaseEncodingGradient(gradSliceRefPhase3D);
        if (isEnablePhase)
            gradPhase2D_comp.refocalizePhaseEncodingGradient(gradPhase2D);
        gradSlicePhase3D_comp.applyAmplitude();
        gradPhase2D_comp.applyAmplitude();
    }

    @Override
    protected void getCrusherGrad() {
        double time_grad_crusher_top = minInstructionDelay;
        double grad_amp_crusher = 0.0;
        double grad_area_crusher = 0.0;

        if (isMultiplanar) {
            grad_area_crusher = getDouble(GRADIENT_AREA_CRUSHER_PI) / ((GradientMath.GAMMA) * sliceThickness);
            getParam(GRADIENT_AREA_CRUSHER).setValue(grad_area_crusher);
            getParam(GRADIENT_CRUSHER_TOP_TIME).setValue(getDouble(GRADIENT_PHASE_APPLICATION_TIME));
            set(Time_grad_crusher_top, GRADIENT_CRUSHER_TOP_TIME);
            grad_amp_crusher = grad_area_crusher / (grad_shape_rise_time + getDouble(GRADIENT_CRUSHER_TOP_TIME)) / gMax * 100.0;
            getParam(GRADIENT_AMP_CRUSHER).setValue(grad_amp_crusher);
            Gradient gradSliceCrusher = Gradient.createGradient(getSequence(), Grad_amp_slice_crush, Time_grad_crusher_top, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
            gradSliceCrusher.addSpoiler(grad_amp_crusher);
            gradSliceCrusher.applyAmplitude();
        } else {
            set(Time_grad_crusher_top, time_grad_crusher_top);
            getParam(GRADIENT_CRUSHER_TOP_TIME).setValue(time_grad_crusher_top);
            getParam(GRADIENT_AMP_CRUSHER).setValue(grad_amp_crusher);
            getParam(GRADIENT_AREA_CRUSHER).setValue(grad_area_crusher);
            getParam(GRADIENT_AREA_CRUSHER_PI).setValue(grad_area_crusher * (GradientMath.GAMMA) * sliceThickness);
        }
    }

    @Override
    protected void getSpoilerGrad() {
        // -------------------------------------------------------------------------------------------------
        // Spoiler Gradient
        // -------------------------------------------------------------------------------------------------
        set(Time_grad_spoiler_top, GRADIENT_CRUSHER_END_TOP_TIME);

        Gradient gradSliceSpoiler = Gradient.createGradient(getSequence(), Grad_amp_spoiler_slice, Time_grad_spoiler_top, Grad_shape_rise_up, Grad_shape_rise_down, Time_grad_ramp, nucleus);
        if (getBoolean(GRADIENT_ENABLE_SPOILER)) {
            gradSliceSpoiler.addSpoiler(getDouble(GRADIENT_AMP_SPOILER));
        }
        gradSliceSpoiler.applyAmplitude();
    }

    @Override
    protected void getUPDisp() {
        this.getParam(TX_ATT).setValue(pulseTX.getAttParamValue());            // display PULSE_ATT
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
        // DYNAMIC SEQUENCE
        // Calculate frame acquisition time
        // Calculate delay between 4D acquisition
        //----------------------------------------------------------------------
        Table dyn_delay = setSequenceTableValues(Time_btw_dyn_frames, Order.Four);
        double frame_acquisition_time = nb_scan_1d * nb_scan_3d * nb_scan_2d * tr;
        //int nb_4d_withoutDyn = numberOfInversionRecovery * numberOfInversionRecovery * (isDixon ? 3 : 1) * numberOfTrigger;
        int nb_4d_withoutDyn = nb_scan_4d / nb_dynamic_acquisition;

        double time_4DFramesNoDyn_withoutDynDelay = (frame_acquisition_time * nb_4d_withoutDyn
                + min_flush_delay * (nb_4d_withoutDyn - 1));

        double interval_between_dynFrames_delay = min_flush_delay;
        if (isDynamic) {
            double time_between_dynFrames_min = ceilToSubDecimal(time_4DFramesNoDyn_withoutDynDelay + min_flush_delay, 3);

            //    double time_between_frames = getDouble(DYN_TIME_BTW_FRAMES);
            double time_between_dynFrames = isDynamicMinTime ? time_between_dynFrames_min : getDouble(DYN_TIME_BTW_FRAMES);
            getParam(DYN_TIME_BTW_FRAMES).setValue(time_between_dynFrames_min);

            if (time_between_dynFrames < (time_between_dynFrames_min)) {
                this.notifyOutOfRangeParam(DYN_TIME_BTW_FRAMES, time_between_dynFrames_min, ((NumberParam) getParam(DYN_TIME_BTW_FRAMES)).getMaxValue(), "Minimum frame acquisition time ");
                time_between_dynFrames = time_between_dynFrames_min;
            }

            interval_between_dynFrames_delay = time_between_dynFrames - time_4DFramesNoDyn_withoutDynDelay;
        }
        dyn_delay.add(interval_between_dynFrames_delay);

        // ------------------------------------------------------------------
        // Total Acquisition Time
        // ------------------------------------------------------------------
        double sequenceTime = (time_4DFramesNoDyn_withoutDynDelay + interval_between_dynFrames_delay) * nb_dynamic_acquisition + tr * nb_preScan;
        getParam(SEQUENCE_TIME).setValue(sequenceTime);

    }

    @Override
    protected void getMultiParaList() {
        double frame_acquisition_time = nb_scan_1d * nb_scan_3d * nb_scan_2d * tr;

        int number_of_MultiSeries = 1;
        double time_between_MultiSeries = 0;
        ArrayList<Number> multiseries_valuesList = new ArrayList<>();
        String multiseries_parametername = "";

        if ((!is_FSE_vs_MultiEcho && echoTrainLength != 1)) {
            number_of_MultiSeries = echoTrainLength;
            time_between_MultiSeries = echo_spacing;
            multiseries_parametername = "TE";
            for (int i = 0; i < number_of_MultiSeries; i++) {
                double multiseries_value = roundToDecimal(te + i * te, 5) * 1e3;
                multiseries_valuesList.add(multiseries_value);
            }
        } else if (models.get(InvRec.class).isEnabled() && models.get(InvRec.class).nb_inversionRecovery != 1) {
            number_of_MultiSeries = models.get(InvRec.class).nb_inversionRecovery;
            time_between_MultiSeries = frame_acquisition_time;
            multiseries_parametername = "TI";
            for (int i = 0; i < number_of_MultiSeries; i++) {
                double IR_time = roundToDecimal(models.get(InvRec.class).inversionRecoveryTime.get(i), 5) * 1e3;
                multiseries_valuesList.add(IR_time);
            }
        } else if (models.get(ExtTrig.class).isEnabled() && models.get(ExtTrig.class).nb_trigger != 1) {
            number_of_MultiSeries = models.get(ExtTrig.class).nb_trigger;
            time_between_MultiSeries = frame_acquisition_time;
            multiseries_parametername = "TRIGGER DELAY";
            for (int i = 0; i < number_of_MultiSeries; i++) {
                double multiseries_value = roundToDecimal((models.get(ExtTrig.class)).triggerTime.get(i), 5) * 1e3;
                multiseries_valuesList.add(multiseries_value);
            }
        }
        getParam(MULTISERIES_PARAMETER_VALUE).setValue(multiseries_valuesList);
        getParam(MULTISERIES_PARAMETER_NAME).setValue(multiseries_parametername);

        ArrayList<Number> acquisition_timesList = new ArrayList<>();
        double acqusition_time;
        for (int i = 0; i < nb_dynamic_acquisition; i++) {
            for (int j = 0; j < number_of_MultiSeries; j++) {
                acqusition_time = (i * frame_acquisition_time + j * time_between_MultiSeries);
                if (i > 0) { // only the first dynamic phase has dummy scans
                    acqusition_time = acqusition_time + nb_preScan * tr;
                }
                acqusition_time = roundToDecimal(acqusition_time, 3);
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
        if (grad_rise_time < blankingDelay) {
            double new_grad_rise_time = ceilToSubDecimal(blankingDelay, 5);
            notifyOutOfRangeParam(GRADIENT_RISE_TIME, new_grad_rise_time, ((NumberParam) getParam(GRADIENT_RISE_TIME)).getMaxValue(), "Gradient ramp time too short because of blanking delay ");
            grad_rise_time = new_grad_rise_time;
        }
        grad_rise_time = ceilToGRT(grad_rise_time);

        set(Time_grad_ramp, grad_rise_time);
        getParam(GRADIENT_RISE_TIME).setValue(grad_rise_time);
    }

}
