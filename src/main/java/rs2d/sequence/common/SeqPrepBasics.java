package rs2d.sequence.common;

import rs2d.commons.log.Log;
import rs2d.spinlab.data.transformPlugin.TransformPlugin;
import rs2d.spinlab.instrument.Instrument;
import rs2d.spinlab.instrument.InstrumentTxChannel;
import rs2d.spinlab.sequence.table.Table;
import rs2d.spinlab.sequence.table.Utility;
import rs2d.spinlab.sequenceGenerator.BaseSequenceGenerator;
import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.sequenceGenerator.util.Hardware;
import rs2d.spinlab.tools.param.*;
import rs2d.spinlab.tools.role.RoleEnum;
import rs2d.spinlab.tools.table.Order;
import rs2d.spinlab.tools.utility.Nucleus;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeMap;

/**
 * Abstract Class SeqPrepBasics
 * prep common functions
 * V1.0- 2021-3-16 XG
 *
 * First-level and Second-level Structures
 * are Mandatory for all type of MR sequences
 *
 * Basic Functions/ Methods
 * are Shared for all sequences
 *
 */

public abstract class SeqPrepBasics extends BaseSequenceGenerator {
    protected TreeMap<Double, RFPulse> rfPulses = new TreeMap<>();
    protected TreeMap<String, ModelInterface> models = new TreeMap<>();
    protected List<String> modalNames;
    protected RFPulse pulseTX;
    protected Gradient gradSlice;

    // Constant
    public final static int offset_channel_memory = 512;
    public final static int phase_channel_memory = 512;
    public final static int amp_channel_memory = 2048;
    public final static int loopIndice_memory = 2048;
    public final static double defaultInstructionDelay = 0.000010;     // single instruction minimal duration
    public final static double minInstructionDelay = 0.000005;     // single instruction minimal duration
    public final static int nb_shape_points = 128;

    // Hardware
    protected double blankingDelay;
    protected Nucleus nucleus;
    protected double observeFrequency;
    protected double protonFrequency;
    protected double gMax;
    protected double spectralWidth;

    // Enables
    protected boolean isMultiplanar;
    protected boolean isKSCenterMode;
    protected boolean isEnablePhase;
    protected boolean isEnablePhase3D;
    protected boolean isEnableSlice;
    protected boolean isEnableRead;
    protected boolean isFovDoubled;

    // View
    protected double fov;
    protected double fovPhase;
    protected double fov3d;
    protected double sliceThickness;
    protected double spacingBetweenSlice;

    // Matrix
    protected int acqMatrixDimension1D;
    protected int acqMatrixDimension2D;
    protected int acqMatrixDimension3D;
    protected int acqMatrixDimension4D;
    protected int userMatrixDimension1D;
    protected int userMatrixDimension2D;
    protected int userMatrixDimension3D;
    protected int userMatrixDimension4D;

    // Position
    protected double off_center_distance_1D;
    protected double off_center_distance_2D;
    protected double off_center_distance_3D;

    // Time
    protected double observation_time;
    protected double min_time_per_acq_point;
    protected double tr;
    protected double te;
    protected double echo_spacing;

    // Loops
    protected int nb_scan_1d;
    protected int nb_scan_2d;
    protected int nb_scan_3d;
    protected int nb_scan_4d;
    protected int nb_averages;
    protected int nb_interleaved_slice;
    protected int nb_planar_excitation;
    protected int nb_slices_acquired_in_single_scan;
    protected int nb_shoot_3d;
    protected int nb_preScan;
    protected int nb_satband;
    protected int echoTrainLength;

    // Plugins
    protected TransformPlugin plugin;

    protected SeqPrepBasics() {
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  first-level  structure (general)
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    @Override
    public void init() {}

    public void generate() throws Exception {
        initUserParam();
        this.beforeRouting();
        if (!this.isRouted()) {
            this.route();
            this.initAfterRouting();//init before setup
        }
        //   if (!getBoolean( SETUP_MODE)) {
        this.afterRouting();    //avoid exception during setup
        // }
        this.checkAndFireException();
    }

    protected void initUserParam() {
        isKSCenterMode = getBoolean(CommonUP.KS_CENTER_MODE);
        isEnablePhase3D = !isKSCenterMode && getBoolean(CommonUP.GRADIENT_ENABLE_PHASE_3D);
        isEnablePhase = !isKSCenterMode && getBoolean(CommonUP.GRADIENT_ENABLE_PHASE);
        isEnableSlice = getBoolean(CommonUP.GRADIENT_ENABLE_SLICE);
        isEnableRead = getBoolean(CommonUP.GRADIENT_ENABLE_READ);
        isMultiplanar = getBoolean(CommonUP.MULTI_PLANAR_EXCITATION);
        isFovDoubled = getBoolean(CommonUP.FOV_DOUBLED);

        userMatrixDimension1D = getInt(CommonUP.USER_MATRIX_DIMENSION_1D);
        userMatrixDimension2D = getInt(CommonUP.USER_MATRIX_DIMENSION_2D);
        userMatrixDimension3D = getInt(CommonUP.USER_MATRIX_DIMENSION_3D);
        userMatrixDimension4D = getInt(CommonUP.USER_MATRIX_DIMENSION_4D);
        acqMatrixDimension1D = getInt(CommonUP.ACQUISITION_MATRIX_DIMENSION_1D);
        acqMatrixDimension2D = getInt(CommonUP.ACQUISITION_MATRIX_DIMENSION_2D);
        acqMatrixDimension3D = getInt(CommonUP.ACQUISITION_MATRIX_DIMENSION_3D);
        acqMatrixDimension4D = getInt(CommonUP.ACQUISITION_MATRIX_DIMENSION_4D);

        off_center_distance_1D = getDouble(CommonUP.OFF_CENTER_FIELD_OF_VIEW_1D);
        off_center_distance_2D = getDouble(CommonUP.OFF_CENTER_FIELD_OF_VIEW_2D);
        off_center_distance_3D = getDouble(CommonUP.OFF_CENTER_FIELD_OF_VIEW_3D);

        fov = getDouble(CommonUP.FIELD_OF_VIEW);
        fovPhase = getDouble(CommonUP.FIELD_OF_VIEW_PHASE);
        fov3d = getDouble(CommonUP.FIELD_OF_VIEW_3D);
        sliceThickness = getDouble(CommonUP.SLICE_THICKNESS);
        spacingBetweenSlice = getDouble(CommonUP.SPACING_BETWEEN_SLICE);

        observation_time = getDouble(CommonUP.ACQUISITION_TIME_PER_SCAN);
        tr = getDouble(CommonUP.REPETITION_TIME);
        te = getDouble(CommonUP.ECHO_TIME);
        echo_spacing = getDouble(CommonUP.ECHO_SPACING);

        nb_preScan = getInt(CommonUP.DUMMY_SCAN);
        nb_averages = getInt(CommonUP.NUMBER_OF_AVERAGES);
        nb_shoot_3d = getInt(CommonUP.NUMBER_OF_SHOOT_3D);
        echoTrainLength = getInt(CommonUP.ECHO_TRAIN_LENGTH);

        InstrumentTxChannel txCh = Instrument.instance().getTxChannels().get(getListInt(CommonUP.TX_ROUTE).get(0));
        blankingDelay = Math.max(minInstructionDelay, txCh.getRfAmpChannel().getBlankingDelay());
        getParam(CommonUP.MODALITY).setValue("MRI");
    }

    protected void beforeRouting() throws Exception {
        Log.debug(getClass(), "------------ BEFORE ROUTING -------------");
        // -----------------------------------------------
        // Ini RX parameters : nucleus, RX gain & frequencies
        // -----------------------------------------------
        iniTxRx();

        // -----------------------------------------------
        // Ini Modality class
        // -----------------------------------------------
        iniModel();

        // -----------------------------------------------
        // Ini Imaging related
        // -----------------------------------------------
        iniImaging();
    }

    protected void initAfterRouting() {
        // -----------------------------------------------
        // ini SeqParam basics
        // -----------------------------------------------
        iniSeqParamBasics();

        // -----------------------------------------------
        // ini SeqParam enable gradient lines
        // -----------------------------------------------
        iniSeqParamEnabled();
    }

    protected void afterRouting() throws Exception {
        Log.debug(getClass(), "------------ AFTER ROUTING -------------");
        // -----------------------------------------------
        // prep Imaging related
        // -----------------------------------------------
        prepImaging();

        // -----------------------------------------------
        // prep Models
        // -----------------------------------------------
        prepModels();

        // --------------------------------------------------------------------------------------------------------------------------------------------
        // TIMING --- TIMING --- TIMING --- TIMING --- TIMING --- TIMING --- TIMING --- TIMING --- TIMING --- TIMING --- TIMING --- TIMING --- TIMING
        // --------------------------------------------------------------------------------------------------------------------------------------------
        prepSeqTiming();

        //--------------------------------------------------------------------------------------
        // prep Export DICOM
        //--------------------------------------------------------------------------------------
        prepDicom();

        //--------------------------------------------------------------------------------------
        // prep Comments
        //--------------------------------------------------------------------------------------
        //prepComments();
    }

    protected void checkAndFireException() {
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  second-level  structures
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected void iniTxRx() throws Exception {
    }

    protected void iniModel() throws Exception {
    }

    protected void iniImaging() throws Exception {
    }

    protected void iniSeqParamBasics() {
    }

    protected void iniSeqParamEnabled() {
    }

    protected void prepImaging() throws Exception {
    }

    protected void prepModels() throws Exception {
    }

    protected void prepSeqTiming() throws Exception {
    }

    protected void prepDicom() {
    }

    protected void prepComments() {
    }

    //--------------------------------------------------------------------------------------
    // basic functions
    //--------------------------------------------------------------------------------------
    protected void setModel(List<String> modalNames, SeqPrep seqPrep) throws Exception {
        this.modalNames = modalNames;
        ModelFactory modelFactory = new ModelFactory();

        if (this.modalNames != null) {
            for (String modalName : this.modalNames) {
                ModelInterface eachModel = modelFactory.getModel(modalName, seqPrep);
                eachModel.init();
                models.put(modalName,eachModel);
            }
            for (String modalName : this.modalNames) {
                ModelInterface eachModel = modelFactory.getModel(modalName, seqPrep);
                eachModel.initFinal();
                models.put(modalName,eachModel);
            }
        }
    }

    protected double clcGradEqRiseTime(GeneratorSequenceParamEnum Grad_shape_rise_up, GeneratorSequenceParamEnum Grad_shape_rise_down, double grad_rise_time) {
        double grad_shape_rise_factor_up = Utility.voltageFillingFactor(getSequenceTable(Grad_shape_rise_up));
        double grad_shape_rise_factor_down = Utility.voltageFillingFactor(getSequenceTable(Grad_shape_rise_down));

        return grad_shape_rise_factor_up * grad_rise_time + grad_shape_rise_factor_down * grad_rise_time;        // shape dependant equivalent rise time
    }

    protected void clcFovPhase() {
        fovPhase = (getBoolean(CommonUP.FOV_SQUARE)) ? getDouble(CommonUP.FIELD_OF_VIEW) : getDouble(CommonUP.FIELD_OF_VIEW_PHASE);
        fovPhase = fovPhase > getDouble(CommonUP.FIELD_OF_VIEW) ? getDouble(CommonUP.FIELD_OF_VIEW) : fovPhase;
        getParam(CommonUP.FIELD_OF_VIEW_PHASE).setValue(fovPhase);
        try {
            getParam(CommonUP.PHASE_FIELD_OF_VIEW_RATIO).setValue((fovPhase / getDouble(CommonUP.FIELD_OF_VIEW) * 100.0));    // FOV ratio for display
        } catch (Exception exception) {
            Log.warning(getClass(), "No such user parameter: \"PHASE_FIELD_OF_VIEW_RATIO\"");
        }
        try {
            getParam(CommonUP.FOV_RATIO_PHASE).setValue(Math.round(fovPhase / getDouble(CommonUP.FIELD_OF_VIEW) * 100.0));    // FOV ratio for display
        } catch (Exception exception) {
            Log.warning(getClass(), "No such user parameter: \"FOV_RATIO_PHASE\"");
        }
    }

    protected int floorEven(double value) {
        return (int) Math.floor(Math.round(value) / 2.0) * 2;
    }

    protected void setSquarePixel(boolean square, int dim) {
        if (square) {
            switch (dim) {
                case 2:
                    this.userMatrixDimension2D = (int) Math.round(getInt(CommonUP.USER_MATRIX_DIMENSION_1D) * fovPhase / getDouble(CommonUP.FIELD_OF_VIEW));
                    getParam(CommonUP.USER_MATRIX_DIMENSION_2D).setValue(this.userMatrixDimension2D);
                    getParam(CommonUP.USER_PARTIAL_PHASE).setValue(100.0);
                    break;
                case 3:
                    this.userMatrixDimension3D = (int) Math.round(getInt(CommonUP.USER_MATRIX_DIMENSION_1D) * fov3d / getDouble(CommonUP.FIELD_OF_VIEW));
                    getParam(CommonUP.USER_MATRIX_DIMENSION_3D).setValue(this.userMatrixDimension3D);
                    getParam(CommonUP.USER_PARTIAL_SLICE).setValue(100.0);
                    break;
            }
        }
    }

//    protected int[] satBandPrep(GeneratorParamEnum satbandOrientation, GeneratorParamEnum orientation, GeneratorParamEnum imageOrientationSubject) {
//        int[] position_sli_ph_rea = new int[6];
//
//        boolean cranial = false;
//        boolean caudal = false;
//        boolean anterior = false;
//        boolean posterior = false;
//        boolean right = false;
//        boolean left = false;
//        if ("CRANIAL".equalsIgnoreCase(getText(satbandOrientation))) {
//            cranial = true;
//        } else if ("CAUDAL".equalsIgnoreCase(getText(satbandOrientation))) {
//            caudal = true;
//        } else if ("CRANIAL AND CAUDAL".equalsIgnoreCase(getText(satbandOrientation))) {
//            cranial = true;
//            caudal = true;
//        } else if ("ANTERIOR".equalsIgnoreCase(getText(satbandOrientation))) {
//            anterior = true;
//        } else if ("POSTERIOR".equalsIgnoreCase(getText(satbandOrientation))) {
//            posterior = true;
//        } else if ("ANTERIOR AND POSTERIOR".equalsIgnoreCase(getText(satbandOrientation))) {
//            anterior = true;
//            posterior = true;
//        } else if ("RIGHT".equalsIgnoreCase(getText(satbandOrientation))) {
//            right = true;
//        } else if ("LEFT".equalsIgnoreCase(getText(satbandOrientation))) {
//            left = true;
//        } else if ("RIGHT AND LEFT".equalsIgnoreCase(getText(satbandOrientation))) {
//            right = true;
//            left = true;
//        } else if ("RIGHT AND LEFT".equalsIgnoreCase(getText(satbandOrientation))) {
//            right = true;
//            left = true;
//        } else if ("ALL".equalsIgnoreCase(getText(satbandOrientation))) {
//            cranial = true;
//            caudal = true;
//            anterior = true;
//            posterior = true;
//            right = true;
//            left = true;
//        }
//
//        position_sli_ph_rea[0] = 0;
//        position_sli_ph_rea[1] = 0;
//        position_sli_ph_rea[2] = 0;
//        position_sli_ph_rea[3] = 0;
//        position_sli_ph_rea[4] = 0;
//        position_sli_ph_rea[5] = 0;
//        if ("AXIAL".equalsIgnoreCase(getText(orientation))) {
//            position_sli_ph_rea[0] = cranial ? 1 : 0;
//            position_sli_ph_rea[1] = caudal ? 1 : 0;
//            position_sli_ph_rea[2] = anterior ? 1 : 0;
//            position_sli_ph_rea[3] = posterior ? 1 : 0;
//            position_sli_ph_rea[4] = right ? 1 : 0;
//            position_sli_ph_rea[5] = left ? 1 : 0;
//        } else if ("SAGITTAL".equalsIgnoreCase(getText(orientation))) {
//            position_sli_ph_rea[0] = left ? 1 : 0;
//            position_sli_ph_rea[1] = right ? 1 : 0;
//            position_sli_ph_rea[2] = anterior ? 1 : 0;
//            position_sli_ph_rea[3] = posterior ? 1 : 0;
//            position_sli_ph_rea[4] = cranial ? 1 : 0;
//            position_sli_ph_rea[5] = caudal ? 1 : 0;
//        } else if ("CORONAL".equalsIgnoreCase(getText(orientation))) {
//            position_sli_ph_rea[0] = anterior ? 1 : 0;
//            position_sli_ph_rea[1] = posterior ? 1 : 0;
//            position_sli_ph_rea[2] = right ? 1 : 0;
//            position_sli_ph_rea[3] = left ? 1 : 0;
//            position_sli_ph_rea[4] = cranial ? 1 : 0;
//            position_sli_ph_rea[5] = caudal ? 1 : 0;
//        } else if ("OBLIQUE".equalsIgnoreCase(getText(orientation))) {
//            List<Double> image_orientation = getListDouble(imageOrientationSubject);
//            double[][] dir_ind = new double[3][3];
//            for (int i = 0; i < 3; i++) {
//                dir_ind[0][i] = image_orientation.get(i);
//                dir_ind[1][i] = image_orientation.get(i + 3);
//            }
//            dir_ind[2][0] = dir_ind[0][1] * dir_ind[1][2] - dir_ind[0][2] * dir_ind[1][1];
//            dir_ind[2][1] = dir_ind[0][2] * dir_ind[1][0] - dir_ind[0][0] * dir_ind[1][2];
//            dir_ind[2][2] = dir_ind[0][0] * dir_ind[1][1] - dir_ind[0][1] * dir_ind[1][0];
//            int i, j;
//            int max_index = 0;
//            double norm_vector_re = Math.sqrt(Math.pow(dir_ind[0][0], 2) + Math.pow(dir_ind[0][1], 2) + Math.pow(dir_ind[0][2], 2));
//            double norm_vector_ph = Math.sqrt(Math.pow(dir_ind[1][0], 2) + Math.pow(dir_ind[1][1], 2) + Math.pow(dir_ind[1][2], 2));
//            double norm_vector_sl = Math.sqrt(Math.pow(dir_ind[2][0], 2) + Math.pow(dir_ind[2][1], 2) + Math.pow(dir_ind[2][2], 2));
//            //normalizing vectors
//            dir_ind[0][0] = dir_ind[0][0] / norm_vector_re;
//            dir_ind[0][1] = dir_ind[0][1] / norm_vector_re;
//            dir_ind[0][2] = dir_ind[0][2] / norm_vector_re;
//            dir_ind[1][0] = dir_ind[1][0] / norm_vector_ph;
//            dir_ind[1][1] = dir_ind[1][1] / norm_vector_ph;
//            dir_ind[1][2] = dir_ind[1][2] / norm_vector_ph;
//            dir_ind[2][0] = dir_ind[2][0] / norm_vector_sl;
//            dir_ind[2][1] = dir_ind[2][1] / norm_vector_sl;
//            dir_ind[2][2] = dir_ind[2][2] / norm_vector_sl;
//            for (i = 0; i < 3; i++) {
//                for (j = 0; j < 3; j++) {
//                    System.out.println("dir_ind[" + i + "][" + j + "]" + dir_ind[i][j]);
//                }
//            }
//
//            // System.out.println(" direction index and dir ind:  "+direction_index[2]+" "+dir_ind[0][2]);
//            int[] max_vector = new int[3];
//
//            // read, phase and slice vector which component has the largest value
//            for (i = 0; i < 3; i++) {
//                for (j = 0; j < 3; j++) {
//                    if (Math.abs(dir_ind[i][j]) >= Math.abs(dir_ind[i][max_index])) {
//                        max_index = j;
//                    }
//                }
//                max_vector[i] = max_index; // storing each vector's maximum value index
//                //System.out.println("max_vector["+i+"]"+max_vector[i]);
//            }
//
//            boolean[][] anatomy_to_local_mx = new boolean[6][6];
//
//            for (i = 0; i < 6; i++) {
//                for (j = 0; j < 6; j++) {
//                    anatomy_to_local_mx[i][j] = false;
//                }
//            }
//
//            for (i = 0; i < 3; i++) {
//
//                if (dir_ind[i][max_vector[i]] < 0) {
//                    anatomy_to_local_mx[i][max_vector[i] + 3] = true;
//                    anatomy_to_local_mx[i + 3][max_vector[i]] = true;
//                } else {
//                    anatomy_to_local_mx[i][max_vector[i]] = true;
//                    anatomy_to_local_mx[i + 3][max_vector[i] + 3] = true;
//                }
//            }
//            boolean[] local_vector = new boolean[6];
//
//            local_vector[0] = false;
//            local_vector[1] = false;
//            local_vector[2] = false;
//            local_vector[3] = false;
//            local_vector[4] = false;
//            local_vector[5] = false;
//
//            boolean[] anatomy_vector = new boolean[6];
//
//            anatomy_vector[0] = right;
//            anatomy_vector[1] = posterior;
//            anatomy_vector[2] = caudal;
//            anatomy_vector[3] = left;
//            anatomy_vector[4] = anterior;
//            anatomy_vector[5] = cranial;
//
//            boolean sum;
//            for (i = 0; i < 6; i++) {
//                sum = false;
//                for (j = 0; j < 6; j++) {
//                    sum = sum || (anatomy_to_local_mx[i][j] & anatomy_vector[j]);
//                    //	System.out.println("sum= "+sum+" + "+anatomy_to_local_mx[i][j]+"*"+anatomy_vector[j]);
//                }
//                local_vector[i] = sum;
//                // System.out.println("local vector "+local_vector[i]);
//
//            }
//            position_sli_ph_rea[4] = local_vector[0] ? 1 : 0;
//            position_sli_ph_rea[2] = local_vector[1] ? 1 : 0;
//            position_sli_ph_rea[0] = local_vector[2] ? 1 : 0;
//            position_sli_ph_rea[5] = local_vector[3] ? 1 : 0;
//            position_sli_ph_rea[3] = local_vector[4] ? 1 : 0;
//            position_sli_ph_rea[1] = local_vector[5] ? 1 : 0;
//
//            // System.out.println("read+ "+position_sli_ph_rea[4]+" phase+ "+position_sli_ph_rea[2]+" slice+ "+position_sli_ph_rea[0]);
//            // System.out.println("read- "+position_sli_ph_rea[5]+" phase- "+position_sli_ph_rea[3]+" slice- "+position_sli_ph_rea[1]);
//        }
//        boolean is_switch = getBoolean(CommonUP.SWITCH_READ_PHASE);
//        boolean phase_pos_temp = position_sli_ph_rea[2] == 1;
//        boolean phase_neg_temp = position_sli_ph_rea[3] == 1;
//        boolean read_pos_temp = position_sli_ph_rea[4] == 1;
//        boolean read_neg_temp = position_sli_ph_rea[5] == 1;
//        if (is_switch) {
//            position_sli_ph_rea[2] = read_pos_temp ? 1 : 0;
//            position_sli_ph_rea[3] = read_neg_temp ? 1 : 0;
//            position_sli_ph_rea[4] = phase_pos_temp ? 1 : 0;
//            position_sli_ph_rea[5] = phase_neg_temp ? 1 : 0;
//        }
//        return position_sli_ph_rea;
//    }

    public double getTx_bandwidth_factor(GeneratorParamEnum tx_shape, GeneratorParamEnum tx_bandwith_factor_param, GeneratorParamEnum tx_bandwith_factor_param3d) {
//        double tx_bandwidth_factor;
//        String tx_shape_name = pl.getTextParam(tx_shape.name()).getValue();
//
//        List<Double> tx_bandwith_factor_table = pl.getListNumberParam(tx_bandwith_factor_param.name()).getValue().stream().map(Number::doubleValue).collect(Collectors.toList());
//        List<Double> tx_bandwith_factor_3D_table = pl.getListNumberParam(tx_bandwith_factor_param3d.name()).getValue().stream().map(Number::doubleValue).collect(Collectors.toList());

        double tx_bandwidth_factor;
        String tx_shape_name = getText(tx_shape);

        List<Double> tx_bandwith_factor_table = getListDouble(tx_bandwith_factor_param);
        List<Double> tx_bandwith_factor_3D_table = getListDouble(tx_bandwith_factor_param3d);

        if (getBoolean(CommonUP.MULTI_PLANAR_EXCITATION)) {
            if ("GAUSSIAN".equalsIgnoreCase(tx_shape_name)) {
                tx_bandwidth_factor = tx_bandwith_factor_table.get(1);
            } else if ("SINC3".equalsIgnoreCase(tx_shape_name)) {
                tx_bandwidth_factor = tx_bandwith_factor_table.get(2);
            } else if ("SINC5".equalsIgnoreCase(tx_shape_name)) {
                tx_bandwidth_factor = tx_bandwith_factor_table.get(3);
            } else if ("RAMP".equalsIgnoreCase(tx_shape_name)) {
                tx_bandwidth_factor = tx_bandwith_factor_table.get(3);
            } else if ("SLR_8_5152".equalsIgnoreCase(tx_shape_name)) {
                tx_bandwidth_factor = tx_bandwith_factor_table.get(4);
            } else if ("SLR_4_2576".equalsIgnoreCase(tx_shape_name)) {
                tx_bandwidth_factor = tx_bandwith_factor_table.get(5);
            } else {
                tx_bandwidth_factor = tx_bandwith_factor_table.get(0);
            }
        } else {
            if ("GAUSSIAN".equalsIgnoreCase(tx_shape_name)) {
                tx_bandwidth_factor = tx_bandwith_factor_3D_table.get(1);
            } else if ("SINC3".equalsIgnoreCase(tx_shape_name)) {
                tx_bandwidth_factor = tx_bandwith_factor_3D_table.get(2);
            } else if ("SINC5".equalsIgnoreCase(tx_shape_name)) {
                tx_bandwidth_factor = tx_bandwith_factor_3D_table.get(3);
            } else if ("RAMP".equalsIgnoreCase(tx_shape_name)) {
                tx_bandwidth_factor = tx_bandwith_factor_3D_table.get(3);
            } else if ("SLR_8_5152".equalsIgnoreCase(tx_shape_name)) {
                tx_bandwidth_factor = tx_bandwith_factor_3D_table.get(4);
            } else if ("SLR_4_2576".equalsIgnoreCase(tx_shape_name)) {
                tx_bandwidth_factor = tx_bandwith_factor_3D_table.get(5);
            } else {
                tx_bandwidth_factor = tx_bandwith_factor_3D_table.get(0);
            }
        }
        return tx_bandwidth_factor;
    }

    protected double getOff_center_distance_1D_2D_3D(int dim) {
        List<Double> image_orientation = getListDouble(CommonUP.IMAGE_ORIENTATION_SUBJECT);
        double[] direction_index = new double[9];
        direction_index[0] = image_orientation.get(0);
        direction_index[1] = image_orientation.get(1);
        direction_index[2] = image_orientation.get(2);
        direction_index[3] = image_orientation.get(3);
        direction_index[4] = image_orientation.get(4);
        direction_index[5] = image_orientation.get(5);
        direction_index[6] = direction_index[1] * direction_index[5] - direction_index[2] * direction_index[4];
        direction_index[7] = direction_index[2] * direction_index[3] - direction_index[0] * direction_index[5];
        direction_index[8] = direction_index[0] * direction_index[4] - direction_index[1] * direction_index[3];

        double norm_vector_read = Math.sqrt(Math.pow(direction_index[0], 2) + Math.pow(direction_index[1], 2) + Math.pow(direction_index[2], 2));
        double norm_vector_phase = Math.sqrt(Math.pow(direction_index[3], 2) + Math.pow(direction_index[4], 2) + Math.pow(direction_index[5], 2));
        double norm_vector_slice = Math.sqrt(Math.pow(direction_index[6], 2) + Math.pow(direction_index[7], 2) + Math.pow(direction_index[8], 2));

        //Offset according to animal position
        double off_center_distance_Z = getDouble(CommonUP.OFF_CENTER_FIELD_OF_VIEW_Z);
        double off_center_distance_Y = getDouble(CommonUP.OFF_CENTER_FIELD_OF_VIEW_Y);
        double off_center_distance_X = getDouble(CommonUP.OFF_CENTER_FIELD_OF_VIEW_X);

        //Offset according to READ PHASE and SLICE
        double off_center_distance;
        switch (dim) {
            case 1:
                off_center_distance = off_center_distance_X * direction_index[0] / norm_vector_read + off_center_distance_Y * direction_index[1] / norm_vector_read + off_center_distance_Z * direction_index[2] / norm_vector_read;
                break;
            case 2:
                off_center_distance = off_center_distance_X * direction_index[3] / norm_vector_phase + off_center_distance_Y * direction_index[4] / norm_vector_phase + off_center_distance_Z * direction_index[5] / norm_vector_phase;
                break;
            case 3:
                off_center_distance = off_center_distance_X * direction_index[6] / norm_vector_slice + off_center_distance_Y * direction_index[7] / norm_vector_slice + off_center_distance_Z * direction_index[8] / norm_vector_slice;
                break;
            default:
                off_center_distance = 0;
                break;
        }
        return off_center_distance;
    }

    protected double getOff_center_distance_X_Y_Z(int dim, double off_center_distance_1D,
                                                  double off_center_distance_2D, double off_center_distance_3D) {
        List<Double> image_orientation = getListDouble(CommonUP.IMAGE_ORIENTATION_SUBJECT);
        double[] direction_index = new double[9];
        direction_index[0] = image_orientation.get(0);
        direction_index[1] = image_orientation.get(1);
        direction_index[2] = image_orientation.get(2);
        direction_index[3] = image_orientation.get(3);
        direction_index[4] = image_orientation.get(4);
        direction_index[5] = image_orientation.get(5);
        direction_index[6] = direction_index[1] * direction_index[5] - direction_index[2] * direction_index[4];
        direction_index[7] = direction_index[2] * direction_index[3] - direction_index[0] * direction_index[5];
        direction_index[8] = direction_index[0] * direction_index[4] - direction_index[1] * direction_index[3];

        double norm_vector_read = Math.sqrt(Math.pow(direction_index[0], 2) + Math.pow(direction_index[1], 2) + Math.pow(direction_index[2], 2));
        double norm_vector_phase = Math.sqrt(Math.pow(direction_index[3], 2) + Math.pow(direction_index[4], 2) + Math.pow(direction_index[5], 2));
        double norm_vector_slice = Math.sqrt(Math.pow(direction_index[6], 2) + Math.pow(direction_index[7], 2) + Math.pow(direction_index[8], 2));

        //Offset according to READ PHASE and SLICE
        double off_center_distance;
        switch (dim) {
            case 1:
                off_center_distance = off_center_distance_1D * direction_index[0] / norm_vector_read + off_center_distance_2D * direction_index[3] / norm_vector_phase + off_center_distance_3D * direction_index[6] / norm_vector_slice;
                break;
            case 2:
                off_center_distance = off_center_distance_1D * direction_index[1] / norm_vector_read + off_center_distance_2D * direction_index[4] / norm_vector_phase + off_center_distance_3D * direction_index[7] / norm_vector_slice;
                break;
            case 3:
                off_center_distance = off_center_distance_1D * direction_index[2] / norm_vector_read + off_center_distance_2D * direction_index[5] / norm_vector_phase + off_center_distance_3D * direction_index[8] / norm_vector_slice;
                break;
            default:
                off_center_distance = 0;
                break;
        }
        return off_center_distance;
    }

    /**
     * Find the next inferior integer which can divide the dividend : dividend /
     * -divisor- = integer
     *
     * @param divisor  dividend / DIVISOR = integer
     * @param dividend DIVIDEND / divisor = integer
     * @return Next inferior integer which is a multiple of the dividend
     */
    protected int getInferiorDivisorToGetModulusZero(int divisor, int dividend) {
        boolean exit = true;
        int div;
        int new_divisor;
        do {
            div = (int) Math.ceil(dividend / ((double) divisor));
            new_divisor = (int) Math.floor(dividend / ((double) div));
            if (dividend % new_divisor == 0) {
                exit = false;
            } else {
                divisor = new_divisor;
            }
        } while (exit);
        return new_divisor;
    }

    protected List<RoleEnum> getPluginAccess() {
        return Collections.singletonList(RoleEnum.User);
    }

    /**
     * From 2D x 3D matrix, return the indices of the scaned PE within elliptic: corner not scans
     *
     * @param matrixDimension2D : 2D dimension
     * @param matrixDimension3D : 2D dimension
     * @param isKSCentred       : true : symetric k'space around zero; false: always go trough k0
     * @return matrix with coordinate of the sampled points: [k2D0, k3D0,     k2D1, k3D1,     k2D2, k3D2,     k2D3, k3D3    ..... ]
     */
    protected ArrayList<Integer> trajEllipticTableBuilder(int matrixDimension2D, int matrixDimension3D, boolean isKSCentred) {
        double Center2D, Center3D;
        if (isKSCentred) {
            Center2D = 1 / 2.0;// symetric k'space around zero
            Center3D = 1 / 2.0;// symetric k'space around zero
        } else {
            Center2D = 1 / 2.0 + ((matrixDimension2D + 1) % 2) / (2.0 * ((float) matrixDimension2D - 1));// always go trough k0
            Center3D = 1 / 2.0 + ((matrixDimension3D + 1) % 2) / (2.0 * ((float) matrixDimension3D - 1));// always go trough k0
        }
        Center2D *= matrixDimension2D;
        Center3D *= matrixDimension3D;

        ArrayList<Integer> position = new ArrayList<>();
        for (int j = 0; j < matrixDimension2D; j++) {
            for (int k = 0; k < matrixDimension3D; k++) {
                if (Math.pow((j - Center2D) / Center2D, 2) + Math.pow((k - Center3D) / Center3D, 2) < 1.0) {
                    position.add(j);
                    position.add(k);
                }
            }
        }
        return position;
    }

    /**
     * In Cam3
     * Split the traj table into 2D and 3D scans, according to memory limitation,
     * NTraj + dummy = 2D * 3D
     * will find the 2D, 3D couple that generates the least number of dummies
     * the Dummy are added at the end of traj
     * note: 2048 Gradient table split in 2 x 1024 when updatedim == true
     * <p>
     * In Cam4
     * Memory is increased
     * no split is allowed anymore
     *
     * @param limMin2D min number of nb_2D
     * @param limMax2D max number of nb_2D
     * @param traj     table with indices in
     * @return int[3] = nb_2D ,nb_3D, nd_dummy
     */
    protected int[] getNbScans2D3DForUpdateDimension(int limMin2D, int limMax2D, ArrayList<Integer> traj) throws Exception {
        int[] nb2D_3D_Dummy = new int[3];
        int nbPE = traj.size() / 2;
        if (Hardware.getSequenceCompiler().getVersionID() < 1303) {
            if (nbPE > 2 * limMax2D) {   //  << Todo exact number of gradient left
                int tmp22 = limMax2D;
                int newNbPE = 0;
                for (int i = limMin2D; i < limMax2D; i++) {
                    int tmpAdditionnalDummy = (int) Math.ceil(nbPE * 1f / (i)) * i - nbPE;
                    //  System.out.println(i+"  "+ tmpAdditionnalDummy+"   "+(int) Math.ceil(totalAcquisitionSpoke * 1f / (i)) * i+"   "+(int) Math.ceil(totalAcquisitionSpoke * 1f / (i)) * i);
                    if (tmpAdditionnalDummy <= tmp22) {
                        tmp22 = tmpAdditionnalDummy;
                        newNbPE = (int) Math.ceil(nbPE * 1f / (i)) * i;
                        nb2D_3D_Dummy[0] = i;
                    }
                }
                nb2D_3D_Dummy[2] = newNbPE - nbPE;
                nbPE = newNbPE;
                nb2D_3D_Dummy[1] = newNbPE / nb2D_3D_Dummy[0];
                for (int i = 0; i < nb2D_3D_Dummy[2]; i++) {
                    System.out.println(" Add  additionnalDummy " + nb2D_3D_Dummy[2]);
                    traj.add(traj.get(traj.size() - 1));
                    traj.add(traj.get(traj.size() - 1));
                }
            } else {
                nb2D_3D_Dummy[0] = nbPE;
                nb2D_3D_Dummy[1] = 1;
            }
        } else {
            if (nbPE > 2 * limMax2D) {
                System.out.println(" Error: nbPE is larger than limMax2D in Cam4");
                throw new Exception("ERROR: Memory overflow! nbPE is larger than limMax2D");
            } else {
                nb2D_3D_Dummy[0] = nbPE;
                nb2D_3D_Dummy[1] = 1;
            }
        }
        return nb2D_3D_Dummy;
    }

    //--------------------------------------------------------------------------------------
    // basic methods
    //--------------------------------------------------------------------------------------

    public static double roundToDecimal(double numberToBeRounded, double order) {
        return Math.round(numberToBeRounded * Math.pow(10, order)) / Math.pow(10, order);
    }

    public static double ceilToSubDecimal(double numberToBeRounded, double Order) {
        return Math.ceil(numberToBeRounded * Math.pow(10, Order)) / Math.pow(10, Order);
    }

    protected Table setSequenceTableValues(GeneratorSequenceParamEnum tableName, Order order, double... values) {
        Table table = getSequenceTable(tableName);
        table.clear();
        table.setOrder(order);
        table.setLocked(true);

        for (double value : values) {
            table.add(value);
        }
        return table;
    }
    @Override
    public int getInt(GeneratorParamEnum userParam) {
        return super.getInt(userParam);
    }

    @Override
    public double getDouble(GeneratorParamEnum userParam) {
        return super.getDouble(userParam);
    }

    @Override
    public boolean getBoolean(GeneratorParamEnum userParam) {
        return super.getBoolean(userParam);
    }

    @Override
    public <T extends Param> T getParam(GeneratorParamEnum userParam) {
        return super.getParam(userParam);
    }

    @Override
    public String getText(GeneratorParamEnum userParam) {
        return super.getText(userParam);
    }

    @Override
    public List<Number> getListNumber(GeneratorParamEnum userParam) {
        return super.getListNumber(userParam);
    }

    @Override
    public List<Double> getListDouble(GeneratorParamEnum userParam) {
        return super.getListDouble(userParam);
    }

    @Override
    public List<Integer> getListInt(GeneratorParamEnum userParam) {
        return super.getListInt(userParam);
    }

    @Override
    public List<String> getListText(GeneratorParamEnum userParam) {
        return super.getListText(userParam);
    }

    @Override
    public boolean hasParam(GeneratorParamEnum userParam) {
        return super.hasParam(userParam);
    }

    @Override
    public void set(GeneratorSequenceParamEnum sequenceParam, Object value) {
        super.set(sequenceParam, value);
    }

    @Override
    protected <T extends Param> T getSequenceParam(GeneratorSequenceParamEnum sequenceParam) {
        return super.getSequenceParam(sequenceParam);
    }

    @Override
    public void notifyOutOfRangeParam(GeneratorParamEnum userParam, Number minValue, Number maxValue, String message) {
        super.notifyOutOfRangeParam(userParam, minValue, maxValue, message);
    }
}

