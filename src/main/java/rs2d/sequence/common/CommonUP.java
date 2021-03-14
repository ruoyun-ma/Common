package rs2d.sequence.common;

import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.tools.param.Param;


public enum CommonUP implements GeneratorParamEnum {
    // Hardware related
    NUCLEUS_1,
    OFFSET_FREQ_1,
    BASE_FREQ_1,
    RECEIVER_GAIN,
    RECEIVER_COUNT,
    INTERMEDIATE_FREQUENCY,
    OBSERVED_FREQUENCY,
    OBSERVED_NUCLEUS,
    SPECTRAL_WIDTH,
    SPECTRAL_WIDTH_PER_PIXEL,
    SPECTRAL_WIDTH_OPT,
    TX_ROUTE,

    // Enables
    KS_CENTER_MODE,
    MULTI_PLANAR_EXCITATION,
    GRADIENT_ENABLE_PHASE_3D,
    GRADIENT_ENABLE_PHASE,
    GRADIENT_ENABLE_SLICE,
    GRADIENT_ENABLE_READ,

    // View
    FOV_SQUARE,
    FOV_DOUBLED,
    FIELD_OF_VIEW,
    FIELD_OF_VIEW_PHASE,
    FIELD_OF_VIEW_3D,
    PHASE_FIELD_OF_VIEW_RATIO,
    FOV_RATIO_PHASE,
    SWITCH_READ_PHASE,
    SLICE_THICKNESS,
    SPACING_BETWEEN_SLICE,

    RESOLUTION_FREQUENCY,
    RESOLUTION_PHASE,
    RESOLUTION_SLICE,
    SQUARE_PIXEL,

    // Matrix
    USER_MATRIX_DIMENSION_1D,
    USER_MATRIX_DIMENSION_2D,
    USER_MATRIX_DIMENSION_3D,
    USER_MATRIX_DIMENSION_4D,
    ACQUISITION_MATRIX_DIMENSION_1D,
    ACQUISITION_MATRIX_DIMENSION_2D,
    ACQUISITION_MATRIX_DIMENSION_3D,
    ACQUISITION_MATRIX_DIMENSION_4D,
    USER_PARTIAL_PHASE,
    USER_PARTIAL_SLICE,
    USER_ZERO_FILLING_2D,
    USER_ZERO_FILLING_3D,

    // Position
    OFF_CENTER_FIELD_OF_VIEW_Z,
    OFF_CENTER_FIELD_OF_VIEW_Y,
    OFF_CENTER_FIELD_OF_VIEW_X,
    OFF_CENTER_FIELD_OF_VIEW_1D,
    OFF_CENTER_FIELD_OF_VIEW_2D,
    OFF_CENTER_FIELD_OF_VIEW_3D,
    IMAGE_ORIENTATION_SUBJECT,
    ORIENTATION,

    // Time
    ACQUISITION_TIME_PER_SCAN,
    REPETITION_TIME,
    ECHO_TIME,
    ECHO_SPACING,
    GRADIENT_RISE_TIME,

    // Loops
    DUMMY_SCAN,
    NUMBER_OF_AVERAGES,
    ECHO_TRAIN_LENGTH,
    NUMBER_OF_SHOOT_3D,

    // TX
    TX_SHAPE,
    TX_BANDWIDTH_FACTOR,
    TX_BANDWIDTH_FACTOR_3D,
    TX_LENGTH_90,
    TX_LENGTH,
    FLIP_ANGLE,

    // Miscellaneous
    SEQUENCE_VERSION,
    MODALITY,
    TRANSFORM_PLUGIN,
    TX_AMP_ATT_AUTO;

    @Override
    public Param build() {
        return null;
    }
}
