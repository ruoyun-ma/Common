package model;

import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.Param;
import rs2d.spinlab.tools.table.Order;

import common.*;
import kernel.*;
import static common.CommonUP.*;
import static common.CommonSP.*;

/**
 *  V1.0 - 2021.3 XG
 *
 */

public class FlowComp implements ModelInterface {
    public SeqPrep parent;
    protected static boolean isFlowCompEnabled;

    public Gradient gradReadPrepFlowComp;
    public Gradient gradSliceRefPhase3DFlowComp;
    public Gradient gradPhase2DFlowComp;

    protected enum UP implements GeneratorParamEnum {
        FLOW_COMPENSATION,
        FLOWCOMP_DURATION,
        ;

        @Override
        public Param build() {
            throw new UnsupportedOperationException("Unable to create a Param");
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum {
        Grad_enable_flowcomp,
        Time_grad_ramp_flowcomp,
        Time_grad_top_flowcomp,
        Grad_amp_read_prep_flowcomp,
        Grad_amp_phase_3D_prep_flowcomp,
        Grad_amp_phase_2D_prep_flowcomp,
        ;
    }

    public FlowComp(SeqPrep parent) {
        this.parent = parent;
    }

    @Override
    public void init() {
        isFlowCompEnabled = parent.getBoolean(UP.FLOW_COMPENSATION);
    }

    @Override
    public void initFinal() throws Exception {
        parent.set(SP.Grad_enable_flowcomp, isFlowCompEnabled);
        setSeqParamTime();
        initPulseandGrad();
    }

    @Override
    public void prep() throws Exception {
    }

    @Override
    public void prepFinal() {
        if (isFlowCompEnabled) {
            gradSliceRefPhase3DFlowComp.applyAmplitude(Order.Three);
            gradPhase2DFlowComp.applyAmplitude(Order.Two);
            gradReadPrepFlowComp.applyAmplitude();
        }
        if (isFlowCompEnabled) {
            gradSliceRefPhase3DFlowComp.applyAmplitude(Order.Three);
            gradPhase2DFlowComp.applyAmplitude(Order.Two);
            gradReadPrepFlowComp.applyAmplitude();
        }
    }

    @Override
    public RFPulse getRfPulses() {
        return null;
    }

    @Override
    public double getDuration() {
        double GradDuration = 2 * parent.getSequenceTable(SP.Time_grad_ramp_flowcomp).get(0).doubleValue()
                + parent.getSequenceTable(SP.Time_grad_top_flowcomp).get(0).doubleValue();
        return GradDuration;
    }

    @Override
    public String getName() {
        return "FlowComp";
    }

    @Override
    public boolean isEnabled() {
        return isFlowCompEnabled;
    }

    protected void initPulseandGrad() {
        gradReadPrepFlowComp = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_read_prep_flowcomp,
                SP.Time_grad_top_flowcomp, Grad_shape_rise_up, Grad_shape_rise_down, SP.Time_grad_ramp_flowcomp, parent.nucleus);

        gradSliceRefPhase3DFlowComp = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_phase_3D_prep_flowcomp,
                SP.Time_grad_top_flowcomp, Grad_shape_rise_up, Grad_shape_rise_down, SP.Time_grad_ramp_flowcomp, parent.nucleus);

        gradPhase2DFlowComp = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_phase_2D_prep_flowcomp, SP.Time_grad_top_flowcomp,
                Grad_shape_rise_up, Grad_shape_rise_down, SP.Time_grad_ramp_flowcomp, parent.nucleus);
    }

    protected void setSeqParamTime() {
        parent.set(SP.Time_grad_ramp_flowcomp, isFlowCompEnabled ? parent.getDouble(GRADIENT_RISE_TIME) : parent.minInstructionDelay);
        parent.set(SP.Time_grad_top_flowcomp, isFlowCompEnabled ? parent.getDouble(UP.FLOWCOMP_DURATION) : parent.minInstructionDelay);
    }
}
