package kernel;

import rs2d.commons.log.Log;
import rs2d.spinlab.data.transformPlugin.TransformPlugin;
import rs2d.spinlab.sequence.table.Table;
import rs2d.spinlab.sequence.table.Utility;
import rs2d.spinlab.sequenceGenerator.BaseSequenceGenerator;
import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.sequenceGenerator.util.Hardware;
import rs2d.spinlab.tools.param.Param;
import rs2d.spinlab.tools.param.TextParam;
import rs2d.spinlab.tools.role.RoleEnum;
import rs2d.spinlab.tools.table.Order;
import rs2d.spinlab.tools.utility.Nucleus;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeMap;

import common.*;
import model.*;

import static common.CommonUP.*;
import static common.CommonSP.*;

/**
 * Abstract Class SeqPrepBasics
 * prep common functions
 * V1.0- 2021-3-16 XG
 * <p>
 * First-level and Second-level Structures
 * are Mandatory for all type of MR sequences
 * <p>
 * Basic Functions/ Methods
 * are Shared for all sequences
 */

public abstract class SeqPrepBasics extends BaseSequenceGenerator {
    public final ModelContainer models = new ModelContainer();
    public List<Integer> rfPulsesAtt = new ArrayList<>();
    public List<RFPulse> rfPulses = new ArrayList<>();
    public RFPulse pulseTX;
//    public RFPulse pulseTX90; //TODO:XG: we better unify them in future
    public RFPulse pulseTX180;
    public Gradient gradSlice;
//    public Gradient gradSlice90;
    public Gradient gradSlice180;
    public Gradient gradReadout;

    // Constant
    public final static int offset_channel_memory = 512;
    public final static int phase_channel_memory = 512;
    public final static int amp_channel_memory = 2048;
    public final static int loopIndice_memory = 2048;
    public final static double defaultInstructionDelay = 0.000010;     // single instruction minimal duration
    public final static double minInstructionDelay = 0.000005;     // single instruction minimal duration

    // Hardware
    public double blankingDelay;
    public Nucleus nucleus;
    public double observeFrequency;
    public double protonFrequency;
    public double gMax;
    public double spectralWidth;

    // Enables
    public boolean isMultiplanar;
    public boolean isKSCenterMode;
    public boolean isEnablePhase;
    public boolean isEnablePhase3D;
    public boolean isEnableSlice;
    public boolean isEnableRead;
    public boolean isFovDoubled;

    // View
    public double fov;
    public double fovPhase;
    public double fov3d;
    public double sliceThickness;
    public double spacingBetweenSlice;

    // Matrix
    public int acqMatrixDimension1D;
    public int acqMatrixDimension2D;
    public int acqMatrixDimension3D;
    public int acqMatrixDimension4D;
    public int userMatrixDimension1D;
    public int userMatrixDimension2D;
    public int userMatrixDimension3D;
    public int userMatrixDimension4D;

    // Position
    public double off_center_distance_1D;
    public double off_center_distance_2D;
    public double off_center_distance_3D;

    // Time
    public double observation_time;
    public double min_time_per_acq_point;
    public double tr;
    public double te;
    public double echo_spacing;

    // Loops
    public int nb_scan_1d;
    public int nb_scan_2d;
    public int nb_scan_3d;
    public int nb_scan_4d;
    public int nb_averages;
    public int nb_interleaved_slice;
    public int nb_planar_excitation;
    //public int nb_slices_acquired_in_single_scan;
    public int nb_shoot_3d;
    public int nb_preScan;
    public int nb_echo_4D;
    public int echoTrainLength;

    // Miscellaneous
    public TransformPlugin plugin;
    public List<String> tx_shape;
    public int nb_shape_points;
    public double grad_rise_time;

    public SeqPrepBasics() {
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  first-level  structure (general)
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    @Override
    public void init() {
        super.init();
        iniModels();
    }

    @Override
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

    public abstract void initUserParam();

    public void beforeRouting() throws Exception {
        Log.debug(getClass(), "------------ BEFORE ROUTING -------------");
        // -----------------------------------------------
        // Ini RX parameters : nucleus, RX gain & frequencies
        // -----------------------------------------------
        iniTxRx();

        // -----------------------------------------------
        // Ini Modality class
        // -----------------------------------------------
        iniFinalModels();

        // -----------------------------------------------
        // Ini Imaging related
        // -----------------------------------------------
        iniImaging();
    }

    public void initAfterRouting() {
        // -----------------------------------------------
        // ini SeqParam basics
        // -----------------------------------------------
        iniSeqParamBasics();

        // -----------------------------------------------
        // ini SeqParam enable gradient lines
        // -----------------------------------------------
        iniSeqParamEnabled();
    }

    public void afterRouting() throws Exception {
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

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  second-level  structures
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected abstract void iniModels();

    protected abstract void iniTxRx() throws Exception;

    protected abstract void iniFinalModels() throws Exception;

    protected abstract void iniImaging() throws Exception;

    protected abstract void iniSeqParamBasics();

    protected abstract void iniSeqParamEnabled();

    protected abstract void prepImaging() throws Exception;

    protected abstract void prepModels() throws Exception;

    protected abstract void prepSeqTiming() throws Exception;

    protected abstract void prepDicom();

    protected void prepComments() {
    }

    //--------------------------------------------------------------------------------------
    // basic functions
    //--------------------------------------------------------------------------------------
    public void setModels(List<ModelInterface> modelClass, SeqPrep seqPrep) {
        models.addAll(modelClass);
        for (ModelInterface eachModel : modelClass) {
            eachModel.init(seqPrep);
        }
    }

    public void setSuggestedValFromListString(List<String> tx_shape, boolean bVal, GeneratorParamEnum... UPs) {
        for (GeneratorParamEnum UP : UPs) {
            ((TextParam) getParam(UP)).setSuggestedValues(tx_shape);
            ((TextParam) getParam(UP)).setRestrictedToSuggested(bVal);
        }
    }

    public double clcGradEqRiseTime(GeneratorSequenceParamEnum Grad_shape_rise_up, GeneratorSequenceParamEnum Grad_shape_rise_down, double grad_rise_time) {
        double grad_shape_rise_factor_up = Utility.voltageFillingFactor(getSequenceTable(Grad_shape_rise_up));
        double grad_shape_rise_factor_down = Utility.voltageFillingFactor(getSequenceTable(Grad_shape_rise_down));

        return grad_shape_rise_factor_up * grad_rise_time + grad_shape_rise_factor_down * grad_rise_time;        // shape dependant equivalent rise time
    }

    public void clcFovPhase() {
        fovPhase = (getBoolean(FOV_SQUARE)) ? getDouble(FIELD_OF_VIEW) : getDouble(FIELD_OF_VIEW_PHASE);
        fovPhase = fovPhase > getDouble(FIELD_OF_VIEW) ? getDouble(FIELD_OF_VIEW) : fovPhase;
        getParam(FIELD_OF_VIEW_PHASE).setValue(fovPhase);
        try {
            getParam(PHASE_FIELD_OF_VIEW_RATIO).setValue((fovPhase / getDouble(FIELD_OF_VIEW) * 100.0));    // FOV ratio for display
        } catch (Exception exception) {
            Log.warning(getClass(), "No such user parameter: \"PHASE_FIELD_OF_VIEW_RATIO\"");
        }
        try {
            getParam(FOV_RATIO_PHASE).setValue(Math.round(fovPhase / getDouble(FIELD_OF_VIEW) * 100.0));    // FOV ratio for display
        } catch (Exception exception) {
            Log.warning(getClass(), "No such user parameter: \"FOV_RATIO_PHASE\"");
        }
    }

    public int floorEven(double value) {
        return (int) Math.floor(Math.round(value) / 2.0) * 2;
    }

    public void setSquarePixel(boolean square, int dim) {
        if (square) {
            switch (dim) {
                case 2:
                    this.userMatrixDimension2D = (int) Math.round(getInt(USER_MATRIX_DIMENSION_1D) * fovPhase / getDouble(FIELD_OF_VIEW));
                    getParam(USER_MATRIX_DIMENSION_2D).setValue(this.userMatrixDimension2D);
                    getParam(USER_PARTIAL_PHASE).setValue(100.0);
                    break;
                case 3:
                    this.userMatrixDimension3D = (int) Math.round(getInt(USER_MATRIX_DIMENSION_1D) * fov3d / getDouble(FIELD_OF_VIEW));
                    getParam(USER_MATRIX_DIMENSION_3D).setValue(this.userMatrixDimension3D);
                    getParam(USER_PARTIAL_SLICE).setValue(100.0);
                    break;
            }
        }
    }

    public double getTx_bandwidth_factor(GeneratorParamEnum tx_shape, GeneratorParamEnum tx_bandwith_factor_param, GeneratorParamEnum tx_bandwith_factor_param3d) {
        double tx_bandwidth_factor;
        String tx_shape_name = getText(tx_shape);

        List<Double> tx_bandwith_factor_table = getListDouble(tx_bandwith_factor_param);
        List<Double> tx_bandwith_factor_3D_table = getListDouble(tx_bandwith_factor_param3d);

        if (getBoolean(MULTI_PLANAR_EXCITATION)) {
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

    public double getOff_center_distance_1D_2D_3D(int dim) {
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

    public double getOff_center_distance_X_Y_Z(int dim, double off_center_distance_1D,
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
    public int getInferiorDivisorToGetModulusZero(int divisor, int dividend) {
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
    public ArrayList<Integer> trajEllipticTableBuilder(int matrixDimension2D, int matrixDimension3D, boolean isKSCentred) {
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
    public int[] getNbScans2D3DForUpdateDimension(int limMin2D, int limMax2D, ArrayList<Integer> traj) throws Exception {
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

    public Table setSequenceTableValues(GeneratorSequenceParamEnum tableName, Order order, double... values) {
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
        if (hasParam(userParam))
            return super.getInt(userParam);
        else
            return -1;
    }

    @Override
    public double getDouble(GeneratorParamEnum userParam) {
        if (hasParam(userParam))
            return super.getDouble(userParam);
        else
            return -1.0;
    }

    @Override
    public boolean getBoolean(GeneratorParamEnum userParam) {
        if (hasParam(userParam))
            return super.getBoolean(userParam);
        else
            return false;
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
    public <T extends Param> T getSequenceParam(GeneratorSequenceParamEnum sequenceParam) {
        return super.getSequenceParam(sequenceParam);
    }

    @Override
    public <T extends Table> T getSequenceTable(GeneratorSequenceParamEnum sequenceTable) {
        return super.getSequenceTable(sequenceTable);
    }

    @Override
    public void notifyOutOfRangeParam(GeneratorParamEnum userParam, Number minValue, Number maxValue, String message) {
        super.notifyOutOfRangeParam(userParam, minValue, maxValue, message);
    }
}

