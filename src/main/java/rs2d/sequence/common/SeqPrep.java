package rs2d.sequence.common;

import rs2d.spinlab.sequence.table.Table;
import rs2d.spinlab.sequenceGenerator.BaseSequenceGenerator;
import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.sequenceGenerator.util.Hardware;
import rs2d.spinlab.tools.role.RoleEnum;
import rs2d.spinlab.tools.table.Order;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Abstract Class SeqPrep
 * prep common functions
 * V1.1- 2021-1-11 XG
 *
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

    // U params
    protected GeneratorParamEnum FOV_SQUARE;
    protected GeneratorParamEnum FIELD_OF_VIEW;
    protected GeneratorParamEnum FIELD_OF_VIEW_PHASE;
    protected GeneratorParamEnum PHASE_FIELD_OF_VIEW_RATIO;
    protected GeneratorParamEnum FOV_RATIO_PHASE;
    protected GeneratorParamEnum USER_MATRIX_DIMENSION_1D;
    protected GeneratorParamEnum USER_MATRIX_DIMENSION_2D;
    protected GeneratorParamEnum OFF_CENTER_FIELD_OF_VIEW_Z;
    protected GeneratorParamEnum OFF_CENTER_FIELD_OF_VIEW_Y;
    protected GeneratorParamEnum OFF_CENTER_FIELD_OF_VIEW_X;
    protected GeneratorParamEnum MULTI_PLANAR_EXCITATION;
    protected GeneratorParamEnum IMAGE_ORIENTATION_SUBJECT;
    protected GeneratorParamEnum SWITCH_READ_PHASE;

    public SeqPrep(Class<? extends Enum> userParamClass){
        FOV_SQUARE = (GeneratorParamEnum) Enum.valueOf(userParamClass, "FOV_SQUARE");
        FIELD_OF_VIEW = (GeneratorParamEnum) Enum.valueOf(userParamClass, "FIELD_OF_VIEW");
        FIELD_OF_VIEW_PHASE = (GeneratorParamEnum) Enum.valueOf(userParamClass, "FIELD_OF_VIEW_PHASE");
        PHASE_FIELD_OF_VIEW_RATIO = (GeneratorParamEnum) Enum.valueOf(userParamClass, "PHASE_FIELD_OF_VIEW_RATIO");
        FOV_RATIO_PHASE = (GeneratorParamEnum) Enum.valueOf(userParamClass, "FOV_RATIO_PHASE");
        USER_MATRIX_DIMENSION_1D = (GeneratorParamEnum) Enum.valueOf(userParamClass, "USER_MATRIX_DIMENSION_1D");
        USER_MATRIX_DIMENSION_2D = (GeneratorParamEnum) Enum.valueOf(userParamClass, "USER_MATRIX_DIMENSION_2D");
        OFF_CENTER_FIELD_OF_VIEW_Z = (GeneratorParamEnum) Enum.valueOf(userParamClass, "OFF_CENTER_FIELD_OF_VIEW_Z");
        OFF_CENTER_FIELD_OF_VIEW_Y = (GeneratorParamEnum) Enum.valueOf(userParamClass, "OFF_CENTER_FIELD_OF_VIEW_Y");
        OFF_CENTER_FIELD_OF_VIEW_X = (GeneratorParamEnum) Enum.valueOf(userParamClass, "OFF_CENTER_FIELD_OF_VIEW_X");
        MULTI_PLANAR_EXCITATION = (GeneratorParamEnum) Enum.valueOf(userParamClass, "MULTI_PLANAR_EXCITATION");
        IMAGE_ORIENTATION_SUBJECT = (GeneratorParamEnum) Enum.valueOf(userParamClass, "IMAGE_ORIENTATION_SUBJECT");
        SWITCH_READ_PHASE = (GeneratorParamEnum) Enum.valueOf(userParamClass, "SWITCH_READ_PHASE");
    }


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  general  methodes
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public void prepareFovPhase() {
        fovPhase = (getBoolean(FOV_SQUARE)) ? getDouble(FIELD_OF_VIEW) : getDouble(FIELD_OF_VIEW_PHASE);
        fovPhase = fovPhase > getDouble(FIELD_OF_VIEW) ? getDouble(FIELD_OF_VIEW) : fovPhase;
        getParam(FIELD_OF_VIEW_PHASE).setValue(fovPhase);
        getParam(PHASE_FIELD_OF_VIEW_RATIO).setValue((fovPhase / getDouble(FIELD_OF_VIEW) * 100.0));    // FOV ratio for display
        getParam(FOV_RATIO_PHASE).setValue(Math.round(fovPhase / getDouble(FIELD_OF_VIEW) * 100.0));    // FOV ratio for display
    }

    protected int floorEven(double value) {
        return (int) Math.floor(Math.round(value) / 2.0) * 2;
    }

    protected void setSquarePixel(boolean square) {
        if (square) {
            this.userMatrixDimension2D = (int) Math.round(getInt(USER_MATRIX_DIMENSION_1D) * fovPhase / getDouble(FIELD_OF_VIEW));
            getParam(USER_MATRIX_DIMENSION_2D).setValue(this.userMatrixDimension2D);
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
        boolean is_switch = getBoolean(SWITCH_READ_PHASE);
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

        if (getBoolean(MULTI_PLANAR_EXCITATION)) {
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
        List<Double> image_orientation = getListDouble(IMAGE_ORIENTATION_SUBJECT);
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
        double off_center_distance_Z = getDouble(OFF_CENTER_FIELD_OF_VIEW_Z);
        double off_center_distance_Y = getDouble(OFF_CENTER_FIELD_OF_VIEW_Y);
        double off_center_distance_X = getDouble(OFF_CENTER_FIELD_OF_VIEW_X);

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
        List<Double> image_orientation = getListDouble(IMAGE_ORIENTATION_SUBJECT);
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

    public List<RoleEnum> getPluginAccess() {
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
     *
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
        }else{
            if (nbPE > 2 * limMax2D){
                System.out.println(" Error: nbPE is larger than limMax2D in Cam4");
                throw new Exception("ERROR: Memory overflow! nbPE is larger than limMax2D");
            }else{
                nb2D_3D_Dummy[0] = nbPE;
                nb2D_3D_Dummy[1] = 1;
            }
        }
        return nb2D_3D_Dummy;
    }
}