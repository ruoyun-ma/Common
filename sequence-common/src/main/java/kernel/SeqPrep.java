package kernel;

import rs2d.spinlab.instrument.Instrument;
import rs2d.spinlab.instrument.InstrumentTxChannel;
import rs2d.spinlab.instrument.util.GradientMath;
import rs2d.spinlab.sequenceGenerator.util.GradientRotation;
import rs2d.spinlab.sequenceGenerator.util.Hardware;
import rs2d.spinlab.tools.utility.GradientAxe;
import rs2d.spinlab.tools.utility.Nucleus;

import static java.util.Arrays.asList;

import common.*;
import model.*;

import java.util.Collections;
import java.util.Iterator;
import java.util.Map;

import static common.CommonUP.*;
import static common.CommonSP.*;

/**
 * Abstract Class SeqPrep
 * prep common functions
 * V1.4- 2021-3-16 XG
 * V1.3- 2021-3-10 XG
 * V1.2- 2021-2-26 XG
 * V1.1- 2021-1-11 XG
 * <p>
 * Override some of the Second-level structure
 * Build Third-level and Fourth-level structure
 * and provide it to different MR sequences
 */


public abstract class SeqPrep extends SeqPrepBasics {

    public SeqPrep() {
        tx_shape = asList("HARD", "GAUSSIAN", "SINC3", "SINC5", "SLR_8_5152", "SLR_4_2576");
        nb_shape_points = 128;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                   second-level  structures (override)
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    @Override
    public void initUserParam() {
        isKSCenterMode = getBoolean(KS_CENTER_MODE);
        isEnablePhase3D = !isKSCenterMode && getBoolean(GRADIENT_ENABLE_PHASE_3D);
        isEnablePhase = !isKSCenterMode && getBoolean(GRADIENT_ENABLE_PHASE);
        isEnableSlice = getBoolean(GRADIENT_ENABLE_SLICE);
        isEnableRead = getBoolean(GRADIENT_ENABLE_READ);
        isMultiplanar = getBoolean(MULTI_PLANAR_EXCITATION);
        isFovDoubled = getBoolean(FOV_DOUBLED);

        userMatrixDimension1D = getInt(USER_MATRIX_DIMENSION_1D);
        userMatrixDimension2D = getInt(USER_MATRIX_DIMENSION_2D);
        userMatrixDimension3D = getInt(USER_MATRIX_DIMENSION_3D);
        userMatrixDimension4D = getInt(USER_MATRIX_DIMENSION_4D);
        acqMatrixDimension1D = getInt(ACQUISITION_MATRIX_DIMENSION_1D);
        acqMatrixDimension2D = getInt(ACQUISITION_MATRIX_DIMENSION_2D);
        acqMatrixDimension3D = getInt(ACQUISITION_MATRIX_DIMENSION_3D);
        acqMatrixDimension4D = getInt(ACQUISITION_MATRIX_DIMENSION_4D);
        off_center_distance_1D = getDouble(OFF_CENTER_FIELD_OF_VIEW_1D);
        off_center_distance_2D = getDouble(OFF_CENTER_FIELD_OF_VIEW_2D);
        off_center_distance_3D = getDouble(OFF_CENTER_FIELD_OF_VIEW_3D);

        fov = getDouble(FIELD_OF_VIEW);
        fovPhase = getDouble(FIELD_OF_VIEW_PHASE);
        fov3d = getDouble(FIELD_OF_VIEW_3D);
        sliceThickness = getDouble(SLICE_THICKNESS);
        spacingBetweenSlice = getDouble(SPACING_BETWEEN_SLICE);

        observation_time = getDouble(ACQUISITION_TIME_PER_SCAN);
        tr = getDouble(REPETITION_TIME);
        te = getDouble(ECHO_TIME);
        echo_spacing = getDouble(ECHO_SPACING);
        grad_rise_time = getDouble(GRADIENT_RISE_TIME);

        nb_preScan = getInt(DUMMY_SCAN);
        nb_averages = getInt(NUMBER_OF_AVERAGES);
        nb_shoot_3d = getInt(NUMBER_OF_SHOOT_3D);
        echoTrainLength = hasParam(ECHO_TRAIN_LENGTH)? getInt(ECHO_TRAIN_LENGTH) : 1;

        spectralWidth = getDouble(SPECTRAL_WIDTH);
        InstrumentTxChannel txCh = Instrument.instance().getTxChannels().get(getListInt(TX_ROUTE).get(0));
        blankingDelay = Math.max(minInstructionDelay, txCh.getRfAmpChannel().getBlankingDelay());
        if (hasParam(MODALITY)) {
            getParam(MODALITY).setValue("MRI");
        }
    }

    @Override
    protected void iniModels() {
    }

    @Override
    protected void iniTxRx() throws Exception {
        getParam(MAGNETIC_FIELD_STRENGTH).setDefaultValue(Instrument.instance().getDevices().getMagnet().getField());
        getParam(DIGITAL_FILTER_SHIFT).setDefaultValue(Instrument.instance().getDevices().getCameleon().getAcquDeadPointCount());
        getParam(DIGITAL_FILTER_REMOVED).setDefaultValue(Instrument.instance().getDevices().getCameleon().isRemoveAcquDeadPoint());

        nucleus = Nucleus.getNucleusForName(getText(NUCLEUS_1));
        protonFrequency = Instrument.instance().getDevices().getMagnet().getProtonFrequency();
        observeFrequency = nucleus.getFrequency(protonFrequency) + getDouble(OFFSET_FREQ_1);
        getParam(BASE_FREQ_1).setValue(nucleus.getFrequency(protonFrequency));

        min_time_per_acq_point = Hardware.getSequenceCompiler().getTransfertTimePerDataPt();
        gMax = GradientMath.getMaxGradientStrength();

        set(Rx_gain, RECEIVER_GAIN);
        getParam(RECEIVER_COUNT).setValue(Instrument.instance().getObservableRxs(nucleus).size());

        set(Intermediate_frequency, Instrument.instance().getIfFrequency());
        getParam(INTERMEDIATE_FREQUENCY).setValue(Instrument.instance().getIfFrequency());

        set(Tx_frequency, observeFrequency);
        getParam(OBSERVED_FREQUENCY).setValue(observeFrequency);

        set(Tx_nucleus, nucleus);
        getParam(OBSERVED_NUCLEUS).setValue(nucleus);
    }

    @Override
    protected void iniPreModels() throws Exception {
        for (ModelInterface eachModel : models)
            eachModel.initPre();
    }

    @Override
    protected void iniFinalModels() throws Exception {
        for (ModelInterface eachModel : models) {
            eachModel.initFinal();
            System.out.println(eachModel.getName() + " " + eachModel.isEnabled());
        }
    }

    @Override
    protected void iniImaging() throws Exception {
        // -----------------------------------------------
        // Ini Acquisition Matrix
        // -----------------------------------------------
        iniAcqMat();

        // -----------------------------------------------
        // Ini TransformPlugin
        // -----------------------------------------------
        iniTransformPlugin();

        // -----------------------------------------------
        // Ini Nb XD
        // -----------------------------------------------
        iniScanLoop();

        // -----------------------------------------------
        // Ini SEQ_DESCRIPTION
        // -----------------------------------------------
        iniSeqDisp();

        // -----------------------------------------------
        // Ini Image Orientation
        // -----------------------------------------------
        iniImgOrientation();
    }

    @Override
    protected void prepImaging() throws Exception {
        // -----------------------------------------------
        // prep RF pulse and sliceGrad parameters
        // -----------------------------------------------
        prepRFandSliceGrad();

        // -----------------------------------------------
        // prep Imaging related gradients
        // -----------------------------------------------
        prepImagingGrads();
    }

    @Override
    protected void prepModels() throws Exception {

        for (ModelInterface eachModel : models) {
            for (int icyc = 0; icyc < 2; icyc++) {
                //cyc it 2 times in case of dependencies
                eachModel.prep();
            }
            if (eachModel.getRfPulses() != null) {
                rfPulses.add(eachModel.getRfPulses()); // We need rfPulses because pulses may be overwritten in rfPulsesTree
            }
        }


        // we find MaxPower for all pulses including pulses in both model and imaging parts
        if (getBoolean(TX_AMP_ATT_AUTO)) {
            for (RFPulse eachPulse : rfPulses) {
                if (eachPulse.getFlipAngle() > 0.0)
                    rfPulsesAtt.add(eachPulse.prepAtt(80, getListInt(TX_ROUTE)));
            }
            set(Tx_att, Collections.min(rfPulsesAtt));

            //---- local ATT
            if (!Collections.max(rfPulsesAtt).equals(Collections.min(rfPulsesAtt))) {
                for (RFPulse eachPulse : rfPulses) {
                    if (Math.abs(eachPulse.getAtt() - Collections.max(rfPulsesAtt))
                            < Math.abs(eachPulse.getAtt() - Collections.min(rfPulsesAtt))) {
                        if (eachPulse.getAttOffset() != -1) {
                            eachPulse.setAttOffset(Collections.max(rfPulsesAtt) - Collections.min(rfPulsesAtt));
                        }
                    } else {
                        if (eachPulse.getAttOffset() != -1) {
                            eachPulse.setAttOffset(0);
                        }
                    }
                }
            }

            for (RFPulse eachPulse : rfPulses) {
                eachPulse.prepTxAmp(getListInt(TX_ROUTE));
            }

            getUPDisp();
        }

        for (ModelInterface eachModel : models) {
            eachModel.prepFinal();
        }
    }

    @Override
    protected void prepSeqTiming() throws Exception {
        // ------------------------------------------
        // calculate Time and Delay
        // ------------------------------------------
        getTimeandDelay();

        // ------------------------------------------
        // calculate TR & search for incoherence
        // ------------------------------------------
        getTR();

        // ------------------------------------------------------------------
        // calculate Total Acquisition Time
        // ------------------------------------------------------------------
        getAcqTime();

        // ------------------------------------------------------------------
        // calculate Acquisition Time Offset and MultiSeries Parameter list
        // ------------------------------------------------------------------
        getMultiParaList();
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                 third-level  structures
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected void iniAcqMat() throws Exception {
        getAcq1D();

        // -----------------------------------------------
        // 2nd D managment
        // -----------------------------------------------
        getAcq2D();

        // -----------------------------------------------
        // 3nd D managment
        // -----------------------------------------------
        getAcq3D();

//        // -----------------------------------------------
//        // 3D managment 2/2: MEMORY LIMITATION & Manage PE trajectory
//        // -----------------------------------------------
//        reprepAcq3D();

        // -----------------------------------------------
        // 4D managment:  Dynamic, MultiEcho, External triggering, Multi Echo
        // -----------------------------------------------
        getAcq4D();

        // -----------------------------------------------
        // Image Resolution
        // -----------------------------------------------
        getImgRes();
    }

    protected void iniETLandVFL() throws Exception {
    }

    protected void iniTransformPlugin() throws Exception {
    }

    protected void iniScanLoop() {
    }

    protected void iniSeqDisp() {
    }

    protected void iniImgOrientation() {
        //READ PHASE and SLICE matrix
        off_center_distance_1D = getOff_center_distance_1D_2D_3D(1);
        off_center_distance_2D = getOff_center_distance_1D_2D_3D(2);
        off_center_distance_3D = getOff_center_distance_1D_2D_3D(3);

        //Offset according to ENABLE READ PHASE and SLICE
        off_center_distance_1D = !isEnableRead ? 0 : off_center_distance_1D;
        off_center_distance_2D = !isEnablePhase ? 0 : off_center_distance_2D;

        if (!isEnableSlice && isMultiplanar || !isEnablePhase3D) {
            off_center_distance_3D = 0;
        }

        getParam(OFF_CENTER_FIELD_OF_VIEW_X).setValue(roundToDecimal(getOff_center_distance_X_Y_Z(1, off_center_distance_1D, off_center_distance_2D, off_center_distance_3D), 5));
        getParam(OFF_CENTER_FIELD_OF_VIEW_Y).setValue(roundToDecimal(getOff_center_distance_X_Y_Z(2, off_center_distance_1D, off_center_distance_2D, off_center_distance_3D), 5));
        getParam(OFF_CENTER_FIELD_OF_VIEW_Z).setValue(roundToDecimal(getOff_center_distance_X_Y_Z(3, off_center_distance_1D, off_center_distance_2D, off_center_distance_3D), 5));

        boolean is_read_phase_inverted = getBoolean(SWITCH_READ_PHASE);
        if (is_read_phase_inverted) {
            set(Gradient_axe_phase, GradientAxe.R);
            set(Gradient_axe_read, GradientAxe.P);
            double off_center_distance_tmp = off_center_distance_2D;
            off_center_distance_2D = off_center_distance_1D;
            off_center_distance_1D = off_center_distance_tmp;
        } else {
            set(Gradient_axe_phase, GradientAxe.P);
            set(Gradient_axe_read, GradientAxe.R);
        }

//        //XG:for 3D PE
//        if (hasParam(SWITCH_READ_SLICE)) {
//            boolean is_read_slice_inverted = getBoolean(SWITCH_READ_SLICE);
//            if (is_read_slice_inverted) {
//                set(Gradient_axe_slice, GradientAxe.R);
//                set(Gradient_axe_read, GradientAxe.S);
//                double off_center_distance_tmp = off_center_distance_3D;
//                off_center_distance_3D = off_center_distance_1D;
//                off_center_distance_1D = off_center_distance_tmp;
//            } else {
//                set(Gradient_axe_slice, GradientAxe.S);
//                set(Gradient_axe_read, GradientAxe.R);
//            }
//        }

        getParam(OFF_CENTER_FIELD_OF_VIEW_3D).setValue(off_center_distance_3D);
        getParam(OFF_CENTER_FIELD_OF_VIEW_2D).setValue(off_center_distance_2D);
        getParam(OFF_CENTER_FIELD_OF_VIEW_1D).setValue(off_center_distance_1D);

        // -----------------------------------------------
        // activate gradient rotation matrix
        // -----------------------------------------------
        GradientRotation.setSequenceGradientRotation(this);
    }

    protected void prepRFandSliceGrad() throws Exception {
    }

    protected void prepImagingGrads() throws Exception {
        // -----------------------------------------------
        // calculate ADC observation time
        // -----------------------------------------------
        getADC();

        // -----------------------------------------------
        // calculate READ gradient amplitude
        // -----------------------------------------------
        getROGrad();

        //----------------------------------------------------------------------
        // OFF CENTER FIELD OF VIEW 1D
        // modify RX FREQUENCY OFFSET
        //----------------------------------------------------------------------
        getRx();

        //--------------------------------------------------------------------------------------
        // PHASE CYCLING
        //--------------------------------------------------------------------------------------
        getPhaseCyc();

        // -------------------------------------------------------------------------------------------------
        // calculate READ_PREP  & SLICE_REF
        // -------------------------------------------------------------------------------------------------
        getPrephaseGrad();

        // -------------------------------------------------------------------------------------------------
        // calculate PHASE_3D  & PHASE_2D
        // -------------------------------------------------------------------------------------------------
        getPEGrad();

        // -------------------------------------------------------------------------------------------------
        // calculate time for 2D_PHASE, 3D_PHASE SLICE_REF or READ_PREP
        // -------------------------------------------------------------------------------------------------
        getGradOpt();

        // -------------------------------------------------------------------------------------------------
        // Calculate Crusher Grad AMPLITUDE
        // -------------------------------------------------------------------------------------------------
        getCrusherGrad();

        // -------------------------------------------------------------------------------------------------
        // Spoiler Gradient
        // -------------------------------------------------------------------------------------------------
        getSpoilerGrad();
    }


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                 fourth-level  structures
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected void getAcq1D() throws Exception {
        // FOV
        double fov_eff = isFovDoubled ? (getDouble(FIELD_OF_VIEW) * 2) : getDouble(FIELD_OF_VIEW);

        // Pixel dimension calculation
        acqMatrixDimension1D = getInt(USER_MATRIX_DIMENSION_1D) * (isFovDoubled ? 2 : 1);
        if(hasParam(RESOLUTION_FREQUENCY))
        getParam(RESOLUTION_FREQUENCY).setValue(fov_eff / acqMatrixDimension1D); // frequency true resolution for display

        // MATRIX
        double spectralWidthPerPixel = getDouble(SPECTRAL_WIDTH_PER_PIXEL);
        spectralWidth = isFovDoubled ? (spectralWidth * 2) : spectralWidth;
        spectralWidth = getBoolean(SPECTRAL_WIDTH_OPT) ? spectralWidth : spectralWidthPerPixel * acqMatrixDimension1D;

//        spectralWidth = Hardware.getNearestSpectralWidth(spectralWidth);      // get real spectral width from Chameleon
        spectralWidth = Hardware.getSequenceCompiler().getNearestSW(spectralWidth);      //  to be replaced by above when correctly implemented for Cam4 05/07/2019
        double spectralWidthUP = isFovDoubled ? (spectralWidth / 2) : spectralWidth;
        spectralWidthPerPixel = spectralWidth / acqMatrixDimension1D;
        getParam(SPECTRAL_WIDTH_PER_PIXEL).setValue(spectralWidthPerPixel);
        getParam(SPECTRAL_WIDTH).setValue(spectralWidthUP);
        observation_time = acqMatrixDimension1D / spectralWidth;
        getParam(ACQUISITION_TIME_PER_SCAN).setValue(observation_time);   // display observation time
        getParam(ACQUISITION_MATRIX_DIMENSION_1D).setValue(acqMatrixDimension1D);
    }

    protected void getAcq2D() {
        // FOV
        clcFovPhase();
        setSquarePixel(getBoolean(SQUARE_PIXEL), 2);

        // MATRIX 2D
        if ("Elliptical3D".equalsIgnoreCase((String) (getParam(TRANSFORM_PLUGIN).getValue()))) {
            if (getDouble(USER_PARTIAL_PHASE) <= 55)
                getParam(USER_PARTIAL_PHASE).setValue(55.0); // for elliptical3D, minimum of 2D partialFourier set 55%
        }

        double partial_phase = getDouble(USER_PARTIAL_PHASE);
        acqMatrixDimension2D = floorEven(partial_phase / 100f * userMatrixDimension2D);
        acqMatrixDimension2D = (acqMatrixDimension2D < 4) && isEnablePhase ? 4 : acqMatrixDimension2D;

        getParam(USER_ZERO_FILLING_2D).setValue((100 - partial_phase));
        getParam(ACQUISITION_MATRIX_DIMENSION_2D).setValue(acqMatrixDimension2D);
    }

    protected void getAcq3D() {
    }

    protected void regetAcq3D() {
    }

    protected void getAcq4D() {
    }

    protected void getImgRes() {
        // PIXEL dimension calculation
        double pixelDimensionPhase = fovPhase / getInt(ACQUISITION_MATRIX_DIMENSION_2D);
        getParam(RESOLUTION_PHASE).setValue(pixelDimensionPhase); // phase true resolution for display

        double pixel_dimension_3D;
        if (isMultiplanar) {
            pixel_dimension_3D = getDouble(SLICE_THICKNESS);
        } else {
            pixel_dimension_3D = getDouble(SLICE_THICKNESS) * getInt(USER_MATRIX_DIMENSION_3D) / getInt(ACQUISITION_MATRIX_DIMENSION_3D); //true resolution
        }
        getParam(RESOLUTION_SLICE).setValue(pixel_dimension_3D); // phase true resolution for display
    }

    protected void getADC() {
        set(Time_rx, getDouble(ACQUISITION_TIME_PER_SCAN));
        set(LO_att, Instrument.instance().getLoAttenuation());
    }

    protected void getROGrad() throws Exception {
    }

    protected void getRx() {
    }

    protected void getPhaseCyc() {
    }

    protected void getPrephaseGrad() {
    }

    protected void getPEGrad() {
    }

    protected void getGradOpt() {
    }

    protected void getCrusherGrad() {
    }

    protected void getSpoilerGrad() {
    }

    //    protected void getFlybackGrad() {
//    }
//
//    protected void getExtTrig() throws Exception {
//        if (this.modalNames != null && modalNames.contains("ExtTrig")) {
//            models.get("ExtTrig").prep();
//        }
//    }
//
//    protected void getIR() {
//    }
//
//    protected void getFatSatGrad() {
//    }
//
//    protected void getFlow() {
//    }
//
//    protected void getSatBandGrad() {
//    }
//
//    protected void getModelsFinal() throws Exception {
//        if (this.modalNames != null && modalNames.contains("ExtTrig")) {
//            models.get("ExtTrig").prepFinal();
//        }
//    }
    protected void getUPDisp() {
    }

    protected void getTimeandDelay() throws Exception {
    }

    protected void getTR() {
    }

    protected void getAcqTime() {
    }

    protected void getMultiParaList() {
    }

}