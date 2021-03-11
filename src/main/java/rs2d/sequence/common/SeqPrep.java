package rs2d.sequence.common;

import com.sun.java.swing.plaf.windows.WindowsInternalFrameTitlePane;
import com.sun.tools.javac.jvm.Gen;
import rs2d.spinlab.instrument.Instrument;
import rs2d.spinlab.instrument.util.GradientMath;
import rs2d.spinlab.sequence.table.Table;
import rs2d.spinlab.sequence.table.Utility;
import rs2d.spinlab.sequenceGenerator.BaseSequenceGenerator;
import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.sequenceGenerator.util.GradientRotation;
import rs2d.spinlab.sequenceGenerator.util.Hardware;
import rs2d.spinlab.tools.param.MriDefaultParams;
import rs2d.spinlab.tools.param.Param;
import rs2d.spinlab.tools.role.RoleEnum;
import rs2d.spinlab.tools.table.Order;
import rs2d.spinlab.tools.utility.GradientAxe;
import rs2d.spinlab.tools.utility.Nucleus;
import sun.security.provider.PolicyParser;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Abstract Class SeqPrep
 * prep common functions
 * V1.1- 2021-1-11 XG
 * V1.2- 2021-2-26 XG
 * V1.3- 2021-3-10 XG
 * <p>
 * one have to change address of the U & S from the local sequence manually
 * TODO: Fix this dependency
 */
//import static rs2d.sequence.spinecho.U.*;
//import rs2d.sequence.spinecho.S;
//import rs2d.sequence.spinecho.U;

public abstract class SeqPrep extends BaseSequenceGenerator {
    // private params
    private double fovPhase;
    private int userMatrixDimension2D;
    private double off_center_distance_1D;
    private double off_center_distance_2D;
    private double off_center_distance_3D;


    // protected params
    protected Nucleus nucleus;
    protected double protonFrequency;
    protected double gMax;
    protected double observeFrequency;
    protected double min_time_per_acq_point;

    protected boolean isMultiplanar;
    protected boolean isKSCenterMode;
    protected boolean isEnablePhase;
    protected boolean isEnablePhase3D;
    protected boolean isEnableSlice;
    protected boolean isEnableRead;


    // U params
    protected GeneratorParamEnum nucleus_1;
    protected GeneratorParamEnum offset_freq_1;
    protected GeneratorParamEnum base_freq_1;
    protected GeneratorParamEnum receiver_gain;
    protected GeneratorParamEnum receiver_count;
    protected GeneratorParamEnum intermediate_frequency;
    protected GeneratorParamEnum observed_frequency;
    protected GeneratorParamEnum observed_nucleus;
    protected GeneratorParamEnum ks_center_mode;
    protected GeneratorParamEnum gradient_enable_phase_3d;
    protected GeneratorParamEnum gradient_enable_phase;
    protected GeneratorParamEnum gradient_enable_slice;
    protected GeneratorParamEnum gradient_enable_read;
    protected GeneratorParamEnum resolution_phase;
    protected GeneratorParamEnum resolution_slice;

    protected GeneratorParamEnum fov_square;
    protected GeneratorParamEnum field_of_view;
    protected GeneratorParamEnum field_of_view_phase;
    protected GeneratorParamEnum phase_field_of_view_ratio;
    protected GeneratorParamEnum fov_ratio_phase;
    protected GeneratorParamEnum slice_thickness;
    protected GeneratorParamEnum user_matrix_dimension_1d;
    protected GeneratorParamEnum user_matrix_dimension_2d;
    protected GeneratorParamEnum user_matrix_dimension_3d;
    protected GeneratorParamEnum acquisition_matrix_dimension_2d;
    protected GeneratorParamEnum acquisition_matrix_dimension_3d;

    protected GeneratorParamEnum off_center_field_of_view_z;
    protected GeneratorParamEnum off_center_field_of_view_y;
    protected GeneratorParamEnum off_center_field_of_view_x;
    protected GeneratorParamEnum multi_planar_excitation;
    protected GeneratorParamEnum image_orientation_subject;
    protected GeneratorParamEnum switch_read_phase;

    protected SeqPrep(Class<? extends Enum> U) {
        nucleus_1 = initialize(U, "NUCLEUS_1");
        offset_freq_1 = initialize(U, "OFFSET_FREQ_1");
        base_freq_1 = initialize(U,"BASE_FREQ_1");
        receiver_gain = initialize(U, "RECEIVER_GAIN");
        receiver_count = initialize(U, "RECEIVER_COUNT");
        intermediate_frequency = initialize(U,"INTERMEDIATE_FREQUENCY");
        observed_frequency = initialize(U,"OBSERVED_FREQUENCY");
        observed_nucleus = initialize(U,"OBSERVED_NUCLEUS");

        ks_center_mode = initialize(U, "KS_CENTER_MODE");
        gradient_enable_phase_3d = initialize(U, "GRADIENT_ENABLE_PHASE_3D");
        gradient_enable_phase = initialize(U, "GRADIENT_ENABLE_PHASE");
        gradient_enable_slice = initialize(U, "GRADIENT_ENABLE_SLICE");
        gradient_enable_read = initialize(U, "GRADIENT_ENABLE_READ");

        resolution_phase = initialize(U, "RESOLUTION_PHASE");
        resolution_slice = initialize(U, "RESOLUTION_SLICE");

        fov_square = initialize(U, "FOV_SQUARE");
        field_of_view = initialize(U, "FIELD_OF_VIEW");
        field_of_view_phase = initialize(U, "FIELD_OF_VIEW_PHASE");
        phase_field_of_view_ratio = initialize(U, "PHASE_FIELD_OF_VIEW_RATIO");
        fov_ratio_phase = initialize(U, "FOV_RATIO_PHASE");
        slice_thickness = initialize(U,"SLICE_THICKNESS");
        user_matrix_dimension_1d = initialize(U, "USER_MATRIX_DIMENSION_1D");
        user_matrix_dimension_2d = initialize(U, "USER_MATRIX_DIMENSION_2D");
        user_matrix_dimension_3d = initialize(U, "USER_MATRIX_DIMENSION_3D");
        acquisition_matrix_dimension_2d = initialize(U, "ACQUISITION_MATRIX_DIMENSION_2D");
        acquisition_matrix_dimension_3d = initialize(U, "ACQUISITION_MATRIX_DIMENSION_3D");

        off_center_field_of_view_z = initialize(U, "OFF_CENTER_FIELD_OF_VIEW_Z");
        off_center_field_of_view_y = initialize(U, "OFF_CENTER_FIELD_OF_VIEW_Y");
        off_center_field_of_view_x = initialize(U, "OFF_CENTER_FIELD_OF_VIEW_X");
        multi_planar_excitation = initialize(U, "MULTI_PLANAR_EXCITATION");
        image_orientation_subject = initialize(U, "IMAGE_ORIENTATION_SUBJECT");
        switch_read_phase = initialize(U, "SWITCH_READ_PHASE");
    }

    protected GeneratorParamEnum initialize(Class<? extends Enum> userParamClass, String inVar) {
        GeneratorParamEnum outVar;
        try {
            outVar = (GeneratorParamEnum) Enum.valueOf(userParamClass, inVar);
        } catch (IllegalArgumentException exception) {
            outVar = null;
        }
        return outVar;
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  general  structure
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
        isKSCenterMode = getBoolean(ks_center_mode);
        isEnablePhase3D = !isKSCenterMode && getBoolean(gradient_enable_phase_3d);
        isEnablePhase = !isKSCenterMode && getBoolean(gradient_enable_phase);
        isEnableSlice = getBoolean(gradient_enable_slice);
        isEnableRead = getBoolean(gradient_enable_read);
        isMultiplanar = getBoolean(multi_planar_excitation);
    }

    protected void beforeRouting() throws Exception {}

    protected void initAfterRouting() {}

    protected void afterRouting() throws Exception {}

    protected void checkAndFireException(){}

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  general  methods
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected void iniTxRx(GeneratorSequenceParamEnum Rx_gain, GeneratorSequenceParamEnum Intermediate_frequency,
                           GeneratorSequenceParamEnum Tx_frequency, GeneratorSequenceParamEnum Tx_nucleus) throws Exception {
        nucleus = Nucleus.getNucleusForName(getText(nucleus_1));
        protonFrequency = Instrument.instance().getDevices().getMagnet().getProtonFrequency();
        observeFrequency = nucleus.getFrequency(protonFrequency) + getDouble(offset_freq_1);
        getParam(base_freq_1).setValue(nucleus.getFrequency(protonFrequency));

        min_time_per_acq_point = Hardware.getSequenceCompiler().getTransfertTimePerDataPt();
        gMax = GradientMath.getMaxGradientStrength();

        set(Rx_gain, receiver_gain);
        getParam(receiver_count).setValue(Instrument.instance().getObservableRxs(nucleus).size());

        set(Intermediate_frequency, Instrument.instance().getIfFrequency());
        getParam(intermediate_frequency).setValue(Instrument.instance().getIfFrequency());

        set(Tx_frequency, observeFrequency);
        getParam(observed_frequency).setValue(observeFrequency);

        set(Tx_nucleus, nucleus);
        getParam(observed_nucleus).setValue(nucleus);
    }

    protected void iniImgOrientation(GeneratorSequenceParamEnum Gradient_axe_phase, GeneratorSequenceParamEnum Gradient_axe_read, GeneratorParamEnum SWITCH_READ_PHASE) {
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

        getParam(off_center_field_of_view_x).setValue(roundToDecimal(getOff_center_distance_X_Y_Z(1, off_center_distance_1D, off_center_distance_2D, off_center_distance_3D), 5));
        getParam(off_center_field_of_view_y).setValue(roundToDecimal(getOff_center_distance_X_Y_Z(2, off_center_distance_1D, off_center_distance_2D, off_center_distance_3D), 5));
        getParam(off_center_field_of_view_z).setValue(roundToDecimal(getOff_center_distance_X_Y_Z(3, off_center_distance_1D, off_center_distance_2D, off_center_distance_3D), 5));


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
//        getParam(OFF_CENTER_FIELD_OF_VIEW_3D).setValue(off_center_distance_3D);
//        getParam(OFF_CENTER_FIELD_OF_VIEW_2D).setValue(off_center_distance_2D);
//        getParam(OFF_CENTER_FIELD_OF_VIEW_1D).setValue(off_center_distance_1D);

        // -----------------------------------------------
        // activate gradient rotation matrix
        // -----------------------------------------------
        GradientRotation.setSequenceGradientRotation(this);
    }

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

    protected void getAcq1D() throws Exception {}

    protected void getAcq2D() {}

    protected void getAcq3D() {}

    protected void regetAcq3D() {}

    protected void getAcq4D() {}

    protected void getImgRes() {
        // PIXEL dimension calculation
        double pixelDimensionPhase = fovPhase / getInt(acquisition_matrix_dimension_2d);
        getParam(resolution_phase).setValue(pixelDimensionPhase); // phase true resolution for display

        double pixel_dimension_3D;
        if (isMultiplanar) {
            pixel_dimension_3D = getDouble(slice_thickness);
        } else {
            pixel_dimension_3D = getDouble(slice_thickness) * getInt(user_matrix_dimension_3d) / getInt(acquisition_matrix_dimension_3d); //true resolution
        }
        getParam(resolution_slice).setValue(pixel_dimension_3D); // phase true resolution for display
    }

    protected double getGradEqRiseTime(GeneratorSequenceParamEnum Grad_shape_rise_up, GeneratorSequenceParamEnum Grad_shape_rise_down, double grad_rise_time) {
        double grad_shape_rise_factor_up = Utility.voltageFillingFactor(getSequenceTable(Grad_shape_rise_up));
        double grad_shape_rise_factor_down = Utility.voltageFillingFactor(getSequenceTable(Grad_shape_rise_down));

        return grad_shape_rise_factor_up * grad_rise_time + grad_shape_rise_factor_down * grad_rise_time;        // shape dependant equivalent rise time
    }

    protected void getADC(GeneratorSequenceParamEnum Time_rx, GeneratorSequenceParamEnum LO_att, double observation_time) {
        set(Time_rx, observation_time);
        set(LO_att, Instrument.instance().getLoAttenuation());
    }

    protected void prepareFovPhase() {
        fovPhase = (getBoolean(fov_square)) ? getDouble(field_of_view) : getDouble(field_of_view_phase);
        fovPhase = fovPhase > getDouble(field_of_view) ? getDouble(field_of_view) : fovPhase;
        getParam(field_of_view_phase).setValue(fovPhase);
        getParam(phase_field_of_view_ratio).setValue((fovPhase / getDouble(field_of_view) * 100.0));    // FOV ratio for display
        getParam(fov_ratio_phase).setValue(Math.round(fovPhase / getDouble(field_of_view) * 100.0));    // FOV ratio for display
    }

    protected int floorEven(double value) {
        return (int) Math.floor(Math.round(value) / 2.0) * 2;
    }

    protected void setSquarePixel(boolean square) {
        if (square) {
            this.userMatrixDimension2D = (int) Math.round(getInt(user_matrix_dimension_1d) * fovPhase / getDouble(field_of_view));
            getParam(user_matrix_dimension_2d).setValue(this.userMatrixDimension2D);
        }
    }

    protected int[] satBandPrep(GeneratorParamEnum satbandOrientation, GeneratorParamEnum orientation, GeneratorParamEnum imageOrientationSubject) {
        int[] position_sli_ph_rea = new int[6];

        boolean cranial = false;
        boolean caudal = false;
        boolean anterior = false;
        boolean posterior = false;
        boolean right = false;
        boolean left = false;
        if ("CRANIAL".equalsIgnoreCase(getText(satbandOrientation))) {
            cranial = true;
        } else if ("CAUDAL".equalsIgnoreCase(getText(satbandOrientation))) {
            caudal = true;
        } else if ("CRANIAL AND CAUDAL".equalsIgnoreCase(getText(satbandOrientation))) {
            cranial = true;
            caudal = true;
        } else if ("ANTERIOR".equalsIgnoreCase(getText(satbandOrientation))) {
            anterior = true;
        } else if ("POSTERIOR".equalsIgnoreCase(getText(satbandOrientation))) {
            posterior = true;
        } else if ("ANTERIOR AND POSTERIOR".equalsIgnoreCase(getText(satbandOrientation))) {
            anterior = true;
            posterior = true;
        } else if ("RIGHT".equalsIgnoreCase(getText(satbandOrientation))) {
            right = true;
        } else if ("LEFT".equalsIgnoreCase(getText(satbandOrientation))) {
            left = true;
        } else if ("RIGHT AND LEFT".equalsIgnoreCase(getText(satbandOrientation))) {
            right = true;
            left = true;
        } else if ("RIGHT AND LEFT".equalsIgnoreCase(getText(satbandOrientation))) {
            right = true;
            left = true;
        } else if ("ALL".equalsIgnoreCase(getText(satbandOrientation))) {
            cranial = true;
            caudal = true;
            anterior = true;
            posterior = true;
            right = true;
            left = true;
        }

        position_sli_ph_rea[0] = 0;
        position_sli_ph_rea[1] = 0;
        position_sli_ph_rea[2] = 0;
        position_sli_ph_rea[3] = 0;
        position_sli_ph_rea[4] = 0;
        position_sli_ph_rea[5] = 0;
        if ("AXIAL".equalsIgnoreCase(getText(orientation))) {
            position_sli_ph_rea[0] = cranial ? 1 : 0;
            position_sli_ph_rea[1] = caudal ? 1 : 0;
            position_sli_ph_rea[2] = anterior ? 1 : 0;
            position_sli_ph_rea[3] = posterior ? 1 : 0;
            position_sli_ph_rea[4] = right ? 1 : 0;
            position_sli_ph_rea[5] = left ? 1 : 0;
        } else if ("SAGITTAL".equalsIgnoreCase(getText(orientation))) {
            position_sli_ph_rea[0] = left ? 1 : 0;
            position_sli_ph_rea[1] = right ? 1 : 0;
            position_sli_ph_rea[2] = anterior ? 1 : 0;
            position_sli_ph_rea[3] = posterior ? 1 : 0;
            position_sli_ph_rea[4] = cranial ? 1 : 0;
            position_sli_ph_rea[5] = caudal ? 1 : 0;
        } else if ("CORONAL".equalsIgnoreCase(getText(orientation))) {
            position_sli_ph_rea[0] = anterior ? 1 : 0;
            position_sli_ph_rea[1] = posterior ? 1 : 0;
            position_sli_ph_rea[2] = right ? 1 : 0;
            position_sli_ph_rea[3] = left ? 1 : 0;
            position_sli_ph_rea[4] = cranial ? 1 : 0;
            position_sli_ph_rea[5] = caudal ? 1 : 0;
        } else if ("OBLIQUE".equalsIgnoreCase(getText(orientation))) {
            List<Double> image_orientation = getListDouble(imageOrientationSubject);
            double[][] dir_ind = new double[3][3];
            for (int i = 0; i < 3; i++) {
                dir_ind[0][i] = image_orientation.get(i);
                dir_ind[1][i] = image_orientation.get(i + 3);
            }
            dir_ind[2][0] = dir_ind[0][1] * dir_ind[1][2] - dir_ind[0][2] * dir_ind[1][1];
            dir_ind[2][1] = dir_ind[0][2] * dir_ind[1][0] - dir_ind[0][0] * dir_ind[1][2];
            dir_ind[2][2] = dir_ind[0][0] * dir_ind[1][1] - dir_ind[0][1] * dir_ind[1][0];
            int i, j;
            int max_index = 0;
            double norm_vector_re = Math.sqrt(Math.pow(dir_ind[0][0], 2) + Math.pow(dir_ind[0][1], 2) + Math.pow(dir_ind[0][2], 2));
            double norm_vector_ph = Math.sqrt(Math.pow(dir_ind[1][0], 2) + Math.pow(dir_ind[1][1], 2) + Math.pow(dir_ind[1][2], 2));
            double norm_vector_sl = Math.sqrt(Math.pow(dir_ind[2][0], 2) + Math.pow(dir_ind[2][1], 2) + Math.pow(dir_ind[2][2], 2));
            //normalizing vectors
            dir_ind[0][0] = dir_ind[0][0] / norm_vector_re;
            dir_ind[0][1] = dir_ind[0][1] / norm_vector_re;
            dir_ind[0][2] = dir_ind[0][2] / norm_vector_re;
            dir_ind[1][0] = dir_ind[1][0] / norm_vector_ph;
            dir_ind[1][1] = dir_ind[1][1] / norm_vector_ph;
            dir_ind[1][2] = dir_ind[1][2] / norm_vector_ph;
            dir_ind[2][0] = dir_ind[2][0] / norm_vector_sl;
            dir_ind[2][1] = dir_ind[2][1] / norm_vector_sl;
            dir_ind[2][2] = dir_ind[2][2] / norm_vector_sl;
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    System.out.println("dir_ind[" + i + "][" + j + "]" + dir_ind[i][j]);
                }
            }

            // System.out.println(" direction index and dir ind:  "+direction_index[2]+" "+dir_ind[0][2]);
            int[] max_vector = new int[3];

            // read, phase and slice vector which component has the largest value
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    if (Math.abs(dir_ind[i][j]) >= Math.abs(dir_ind[i][max_index])) {
                        max_index = j;
                    }
                }
                max_vector[i] = max_index; // storing each vector's maximum value index
                //System.out.println("max_vector["+i+"]"+max_vector[i]);
            }

            boolean[][] anatomy_to_local_mx = new boolean[6][6];

            for (i = 0; i < 6; i++) {
                for (j = 0; j < 6; j++) {
                    anatomy_to_local_mx[i][j] = false;
                }
            }

            for (i = 0; i < 3; i++) {

                if (dir_ind[i][max_vector[i]] < 0) {
                    anatomy_to_local_mx[i][max_vector[i] + 3] = true;
                    anatomy_to_local_mx[i + 3][max_vector[i]] = true;
                } else {
                    anatomy_to_local_mx[i][max_vector[i]] = true;
                    anatomy_to_local_mx[i + 3][max_vector[i] + 3] = true;
                }
            }
            boolean[] local_vector = new boolean[6];

            local_vector[0] = false;
            local_vector[1] = false;
            local_vector[2] = false;
            local_vector[3] = false;
            local_vector[4] = false;
            local_vector[5] = false;

            boolean[] anatomy_vector = new boolean[6];

            anatomy_vector[0] = right;
            anatomy_vector[1] = posterior;
            anatomy_vector[2] = caudal;
            anatomy_vector[3] = left;
            anatomy_vector[4] = anterior;
            anatomy_vector[5] = cranial;

            boolean sum;
            for (i = 0; i < 6; i++) {
                sum = false;
                for (j = 0; j < 6; j++) {
                    sum = sum || (anatomy_to_local_mx[i][j] & anatomy_vector[j]);
                    //	System.out.println("sum= "+sum+" + "+anatomy_to_local_mx[i][j]+"*"+anatomy_vector[j]);
                }
                local_vector[i] = sum;
                // System.out.println("local vector "+local_vector[i]);

            }
            position_sli_ph_rea[4] = local_vector[0] ? 1 : 0;
            position_sli_ph_rea[2] = local_vector[1] ? 1 : 0;
            position_sli_ph_rea[0] = local_vector[2] ? 1 : 0;
            position_sli_ph_rea[5] = local_vector[3] ? 1 : 0;
            position_sli_ph_rea[3] = local_vector[4] ? 1 : 0;
            position_sli_ph_rea[1] = local_vector[5] ? 1 : 0;

            // System.out.println("read+ "+position_sli_ph_rea[4]+" phase+ "+position_sli_ph_rea[2]+" slice+ "+position_sli_ph_rea[0]);
            // System.out.println("read- "+position_sli_ph_rea[5]+" phase- "+position_sli_ph_rea[3]+" slice- "+position_sli_ph_rea[1]);
        }
        boolean is_switch = getBoolean(switch_read_phase);
        boolean phase_pos_temp = position_sli_ph_rea[2] == 1;
        boolean phase_neg_temp = position_sli_ph_rea[3] == 1;
        boolean read_pos_temp = position_sli_ph_rea[4] == 1;
        boolean read_neg_temp = position_sli_ph_rea[5] == 1;
        if (is_switch) {
            position_sli_ph_rea[2] = read_pos_temp ? 1 : 0;
            position_sli_ph_rea[3] = read_neg_temp ? 1 : 0;
            position_sli_ph_rea[4] = phase_pos_temp ? 1 : 0;
            position_sli_ph_rea[5] = phase_neg_temp ? 1 : 0;
        }
        return position_sli_ph_rea;
    }

    protected double getTx_bandwidth_factor(GeneratorParamEnum tx_shape, GeneratorParamEnum tx_bandwith_factor_param, GeneratorParamEnum tx_bandwith_factor_param3d) {
        double tx_bandwidth_factor;

        List<Double> tx_bandwith_factor_table = getListDouble(tx_bandwith_factor_param);
        List<Double> tx_bandwith_factor_3D_table = getListDouble(tx_bandwith_factor_param3d);

        if (isMultiplanar) {
            if ("GAUSSIAN".equalsIgnoreCase(getText(tx_shape))) {
                tx_bandwidth_factor = tx_bandwith_factor_table.get(1);
            } else if ("SINC3".equalsIgnoreCase(getText(tx_shape))) {
                tx_bandwidth_factor = tx_bandwith_factor_table.get(2);
            } else if ("SINC5".equalsIgnoreCase(getText(tx_shape))) {
                tx_bandwidth_factor = tx_bandwith_factor_table.get(3);
            } else if ("RAMP".equalsIgnoreCase(getText(tx_shape))) {
                tx_bandwidth_factor = tx_bandwith_factor_table.get(3);
            } else if ("SLR_8_5152".equalsIgnoreCase(getText(tx_shape))) {
                tx_bandwidth_factor = tx_bandwith_factor_table.get(4);
            } else if ("SLR_4_2576".equalsIgnoreCase(getText(tx_shape))) {
                tx_bandwidth_factor = tx_bandwith_factor_table.get(5);
            } else {
                tx_bandwidth_factor = tx_bandwith_factor_table.get(0);
            }
        } else {
            if ("GAUSSIAN".equalsIgnoreCase(getText(tx_shape))) {
                tx_bandwidth_factor = tx_bandwith_factor_3D_table.get(1);
            } else if ("SINC3".equalsIgnoreCase(getText(tx_shape))) {
                tx_bandwidth_factor = tx_bandwith_factor_3D_table.get(2);
            } else if ("SINC5".equalsIgnoreCase(getText(tx_shape))) {
                tx_bandwidth_factor = tx_bandwith_factor_3D_table.get(3);
            } else if ("RAMP".equalsIgnoreCase(getText(tx_shape))) {
                tx_bandwidth_factor = tx_bandwith_factor_3D_table.get(3);
            } else if ("SLR_8_5152".equalsIgnoreCase(getText(tx_shape))) {
                tx_bandwidth_factor = tx_bandwith_factor_3D_table.get(4);
            } else if ("SLR_4_2576".equalsIgnoreCase(getText(tx_shape))) {
                tx_bandwidth_factor = tx_bandwith_factor_3D_table.get(5);
            } else {
                tx_bandwidth_factor = tx_bandwith_factor_3D_table.get(0);
            }
        }
        return tx_bandwidth_factor;
    }

    protected double ceilToSubDecimal(double numberToBeRounded, double Order) {
        return Math.ceil(numberToBeRounded * Math.pow(10, Order)) / Math.pow(10, Order);
    }

    protected double roundToDecimal(double numberToBeRounded, double order) {
        return Math.round(numberToBeRounded * Math.pow(10, order)) / Math.pow(10, order);
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

    protected double getOff_center_distance_1D_2D_3D(int dim) {
        List<Double> image_orientation = getListDouble(image_orientation_subject);
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
        double off_center_distance_Z = getDouble(off_center_field_of_view_z);
        double off_center_distance_Y = getDouble(off_center_field_of_view_y);
        double off_center_distance_X = getDouble(off_center_field_of_view_x);

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
        List<Double> image_orientation = getListDouble(image_orientation_subject);
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
}