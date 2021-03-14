package rs2d.sequence.common;

import rs2d.spinlab.instrument.Instrument;
import rs2d.spinlab.instrument.util.GradientMath;
import rs2d.spinlab.sequenceGenerator.util.GradientRotation;
import rs2d.spinlab.sequenceGenerator.util.Hardware;
import rs2d.spinlab.tools.param.NumberParam;
import rs2d.spinlab.tools.param.TextParam;
import rs2d.spinlab.tools.utility.GradientAxe;
import rs2d.spinlab.tools.utility.Nucleus;

import java.util.List;

import static java.util.Arrays.asList;

/**
 * Abstract Class SeqPrep
 * prep common functions
 * V1.4- 2021-3-16 XG
 * V1.3- 2021-3-10 XG
 * V1.2- 2021-2-26 XG
 * V1.1- 2021-1-11 XG
 *
 * Override some of the Second-level structure
 * Build Third-level and Fourth-level structure
 * and provide it to different MR sequences
 *
 */


public abstract class SeqPrep extends SeqPrepBasics {
    protected SeqPrep() {}

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                   second-level  structures (override)
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    @Override
    public void init() {
        super.init();
        List<String> tx_shape = asList(
                "HARD",
                "GAUSSIAN",
                "SINC3",
                "SINC5",
                "SLR_8_5152",
                "SLR_4_2576");
        ((TextParam) getParam(CommonUP.TX_SHAPE)).setSuggestedValues(tx_shape);
        ((TextParam) getParam(CommonUP.TX_SHAPE)).setRestrictedToSuggested(true);
    }

    @Override
    protected void iniTxRx() throws Exception {
        nucleus = Nucleus.getNucleusForName(getText(CommonUP.NUCLEUS_1));
        protonFrequency = Instrument.instance().getDevices().getMagnet().getProtonFrequency();
        observeFrequency = nucleus.getFrequency(protonFrequency) + getDouble(CommonUP.OFFSET_FREQ_1);
        getParam(CommonUP.BASE_FREQ_1).setValue(nucleus.getFrequency(protonFrequency));

        min_time_per_acq_point = Hardware.getSequenceCompiler().getTransfertTimePerDataPt();
        gMax = GradientMath.getMaxGradientStrength();

        set(CommonSP.Rx_gain, CommonUP.RECEIVER_GAIN);
        getParam(CommonUP.RECEIVER_COUNT).setValue(Instrument.instance().getObservableRxs(nucleus).size());

        set(CommonSP.Intermediate_frequency, Instrument.instance().getIfFrequency());
        getParam(CommonUP.INTERMEDIATE_FREQUENCY).setValue(Instrument.instance().getIfFrequency());

        set(CommonSP.Tx_frequency, observeFrequency);
        getParam(CommonUP.OBSERVED_FREQUENCY).setValue(observeFrequency);

        set(CommonSP.Tx_nucleus, nucleus);
        getParam(CommonUP.OBSERVED_NUCLEUS).setValue(nucleus);
    }

    @Override
    protected void iniImaging() throws Exception {
        // -----------------------------------------------
        // Ini Acquisition Matrix
        // -----------------------------------------------
        iniAcqMat();

        // -----------------------------------------------
        // Ini ETL and VFL
        // -----------------------------------------------
        iniETLandVFL();

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

//    @Override
//    protected void prepModels() throws Exception {
//        //--------------------------------------------------------------------------------------
//        // get Flyback gradient
//        //--------------------------------------------------------------------------------------
//        getFlybackGrad();
//
//        //--------------------------------------------------------------------------------------
//        // get External triggering
//        //--------------------------------------------------------------------------------------
//        getExtTrig();
//
//        // -------------------------------------------------------------------------------------------------
//        // get Inversion Recovery gradients
//        // -------------------------------------------------------------------------------------------------
//        getIR();
//
//        //--------------------------------------------------------------------------------------
//        // get Fat Sat gradients
//        //--------------------------------------------------------------------------------------
//        getFatSatGrad();
//
//        // ------------------------------------------------------------------
//        // get blanking smartTTL_FatSat_table
//        // ------------------------------------------------------------------
//        //prepFatSatSmartTTL();
//
//        //--------------------------------------------------------------------------------------
//        // get Flow
//        //--------------------------------------------------------------------------------------
//        getFlow();
//
//        //--------------------------------------------------------------------------------------
//        // get Sat Band gradients
//        //--------------------------------------------------------------------------------------
//        getSatBandGrad();
//
//        //--------------------------------------------------------------------------------------
//        // get Models Final
//        //--------------------------------------------------------------------------------------
//        getModelsFinal();
//    }

    @Override
    protected void prepModels() throws Exception {
        if (modalNames != null) {
            for (String modalName : modalNames) {
                models.get(modalName).prep();
                if (models.get(modalName).isEnabled()) {
                    if (models.get(modalName).getRfPulses() != null)
                        rfPulses.put(models.get(modalName).getRfPulses().getPower(), models.get(modalName).getRfPulses());
                }
            }

            if (getBoolean(CommonUP.TX_AMP_ATT_AUTO)) {
                RFPulse pulseMaxPower = rfPulses.get(rfPulses.lastKey());
                pulseMaxPower.prepAtt(80, getListInt(CommonUP.TX_ROUTE));

                for (Double key : rfPulses.keySet()) {
                    System.out.println("key "+key + " power "+ rfPulses.get(key).toString());

                    rfPulses.get(key).prepTxAmp(getListInt(CommonUP.TX_ROUTE));
                }
            }

            for (String modalName : modalNames) {
                models.get(modalName).prepFinal();
            }

            getUPDisp();
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

        getParam(CommonUP.OFF_CENTER_FIELD_OF_VIEW_X).setValue(roundToDecimal(getOff_center_distance_X_Y_Z(1, off_center_distance_1D, off_center_distance_2D, off_center_distance_3D), 5));
        getParam(CommonUP.OFF_CENTER_FIELD_OF_VIEW_Y).setValue(roundToDecimal(getOff_center_distance_X_Y_Z(2, off_center_distance_1D, off_center_distance_2D, off_center_distance_3D), 5));
        getParam(CommonUP.OFF_CENTER_FIELD_OF_VIEW_Z).setValue(roundToDecimal(getOff_center_distance_X_Y_Z(3, off_center_distance_1D, off_center_distance_2D, off_center_distance_3D), 5));


        boolean is_read_phase_inverted = getBoolean(CommonUP.SWITCH_READ_PHASE);
        if (is_read_phase_inverted) {
            set(CommonSP.Gradient_axe_phase, GradientAxe.R);
            set(CommonSP.Gradient_axe_read, GradientAxe.P);
            double off_center_distance_tmp = off_center_distance_2D;
            off_center_distance_2D = off_center_distance_1D;
            off_center_distance_1D = off_center_distance_tmp;
        } else {
            set(CommonSP.Gradient_axe_phase, GradientAxe.P);
            set(CommonSP.Gradient_axe_read, GradientAxe.R);
        }
        getParam(CommonUP.OFF_CENTER_FIELD_OF_VIEW_3D).setValue(off_center_distance_3D);
        getParam(CommonUP.OFF_CENTER_FIELD_OF_VIEW_2D).setValue(off_center_distance_2D);
        getParam(CommonUP.OFF_CENTER_FIELD_OF_VIEW_1D).setValue(off_center_distance_1D);

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
        double fov_eff = isFovDoubled ? (getDouble(CommonUP.FIELD_OF_VIEW) * 2) : getDouble(CommonUP.FIELD_OF_VIEW);

        // Pixel dimension calculation
        acqMatrixDimension1D = getInt(CommonUP.USER_MATRIX_DIMENSION_1D) * (isFovDoubled ? 2 : 1);
        getParam(CommonUP.RESOLUTION_FREQUENCY).setValue(fov_eff / acqMatrixDimension1D); // frequency true resolution for display

        // MATRIX
        double spectralWidthPerPixel = getDouble(CommonUP.SPECTRAL_WIDTH_PER_PIXEL);
        spectralWidth = isFovDoubled ? (spectralWidth * 2) : spectralWidth;
        spectralWidth = getBoolean(CommonUP.SPECTRAL_WIDTH_OPT) ? spectralWidth : spectralWidthPerPixel * acqMatrixDimension1D;

//        spectralWidth = Hardware.getNearestSpectralWidth(spectralWidth);      // get real spectral width from Chameleon
        spectralWidth = Hardware.getSequenceCompiler().getNearestSW(spectralWidth);      //  to be replaced by above when correctly implemented for Cam4 05/07/2019
        double spectralWidthUP = isFovDoubled ? (spectralWidth / 2) : spectralWidth;
        spectralWidthPerPixel = spectralWidth / acqMatrixDimension1D;
        getParam(CommonUP.SPECTRAL_WIDTH_PER_PIXEL).setValue(spectralWidthPerPixel);
        getParam(CommonUP.SPECTRAL_WIDTH).setValue(spectralWidthUP);
        observation_time = acqMatrixDimension1D / spectralWidth;
        getParam(CommonUP.ACQUISITION_TIME_PER_SCAN).setValue(observation_time);   // display observation time
        getParam(CommonUP.ACQUISITION_MATRIX_DIMENSION_1D).setValue(acqMatrixDimension1D);
    }

    protected void getAcq2D() {
        // FOV
        clcFovPhase();
        setSquarePixel(getBoolean(CommonUP.SQUARE_PIXEL), 2);

        // MATRIX 2D
        if ("Elliptical3D".equalsIgnoreCase((String) (getParam(CommonUP.TRANSFORM_PLUGIN).getValue()))) {
            if (getDouble(CommonUP.USER_PARTIAL_PHASE) <= 55)
                getParam(CommonUP.USER_PARTIAL_PHASE).setValue(55.0); // for elliptical3D, minimum of 2D partialFourier set 55%
        }
        double partial_phase = getDouble(CommonUP.USER_PARTIAL_PHASE);
        acqMatrixDimension2D = floorEven(partial_phase / 100f * userMatrixDimension2D);
        acqMatrixDimension2D = (acqMatrixDimension2D < 4) && isEnablePhase ? 4 : acqMatrixDimension2D;
        getParam(CommonUP.USER_ZERO_FILLING_2D).setValue((100 - partial_phase));
        getParam(CommonUP.ACQUISITION_MATRIX_DIMENSION_2D).setValue(acqMatrixDimension2D);
    }

    protected void getAcq3D() {
        userParams.addParam(new NumberParam());
    }

    protected void regetAcq3D() {
    }

    protected void getAcq4D() {
    }

    protected void getImgRes() {
        // PIXEL dimension calculation
        double pixelDimensionPhase = fovPhase / getInt(CommonUP.ACQUISITION_MATRIX_DIMENSION_2D);
        getParam(CommonUP.RESOLUTION_PHASE).setValue(pixelDimensionPhase); // phase true resolution for display

        double pixel_dimension_3D;
        if (isMultiplanar) {
            pixel_dimension_3D = getDouble(CommonUP.SLICE_THICKNESS);
        } else {
            pixel_dimension_3D = getDouble(CommonUP.SLICE_THICKNESS) * getInt(CommonUP.USER_MATRIX_DIMENSION_3D) / getInt(CommonUP.ACQUISITION_MATRIX_DIMENSION_3D); //true resolution
        }
        getParam(CommonUP.RESOLUTION_SLICE).setValue(pixel_dimension_3D); // phase true resolution for display
    }

    protected void getADC() {
        set(CommonSP.Time_rx, getDouble(CommonUP.ACQUISITION_TIME_PER_SCAN));
        set(CommonSP.LO_att, Instrument.instance().getLoAttenuation());
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