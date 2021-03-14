package rs2d.sequence.common;

import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.Param;

public enum CommonSP implements GeneratorSequenceParamEnum{
    Tx_att,
    Rx_gain,
    Intermediate_frequency,
    Tx_frequency,
    Tx_nucleus,
    Time_rx,
    LO_att,
    Grad_shape_rise_up,
    Grad_shape_rise_down,
    Gradient_axe_phase,
    Gradient_axe_read;
}
