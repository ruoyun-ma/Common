package model;

import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.Param;

import common.*;
import kernel.*;
import static common.CommonUP.*;
import static common.CommonSP.*;

/**
 *  V1.0 - 2021.3 XG
 *  We decide not to make FlyBack as an model for now
 *
 */

public class FlyBack implements ModelInterface {
    public SeqPrep parent;
    protected static boolean isFlyBackEnabled;

    public Gradient gradReadoutFlyback;
    public double timeGradTopFlyback;
    public double time_flyback_ramp;

    protected enum UP implements GeneratorParamEnum {
        FLYBACK,
        GRADIENT_FLYBACK_TIME,
        ;

        @Override
        public Param build() {
            throw new UnsupportedOperationException("Unable to create a Param");
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum {
        Grad_enable_flyback,
        Grad_amp_flyback,
        Time_flyback,
        Time_grad_ramp_flyback,
        ;
    }

    public FlyBack() {
    }

    @Override
    public void init(SeqPrep parent) {
        this.parent = parent;
    }

    @Override
    public void initPre() throws Exception {
        isFlyBackEnabled = parent.getBoolean(UP.FLYBACK);
    }

    @Override
    public void initFinal() throws Exception {
        parent.set(SP.Grad_enable_flyback, isFlyBackEnabled);

        setSeqParamTime();
        initPulseandGrad();
    }

    @Override
    public void prep() throws Exception {
        prepGrad();
    }

    @Override
    public void prepFinal() {

    }

    @Override
    public RFPulse getRfPulses() {
        return null;
    }

    @Override
    public double getDuration() {
        double GradDuration = 2 * parent.getSequenceTable(SP.Time_grad_ramp_flyback).get(0).doubleValue()
                + parent.getSequenceTable(SP.Time_flyback).get(0).doubleValue();

        return GradDuration;
    }

    @Override
    public String getName() {
        return "FlyBack";
    }

    @Override
    public boolean isEnabled() {
        return isFlyBackEnabled;
    }

    protected void initPulseandGrad() throws Exception {
        gradReadoutFlyback = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_flyback, SP.Time_flyback,
                Grad_shape_rise_up, Grad_shape_rise_down, SP.Time_grad_ramp_flyback, parent.nucleus);
    }

    protected void prepGrad() {
        if (isFlyBackEnabled) {
            double grad_shape_rise_time = parent.clcGradEqRiseTime(Grad_shape_rise_up, Grad_shape_rise_down, parent.getDouble(GRADIENT_RISE_TIME));
            gradReadoutFlyback.refocalizeTotalGradient(parent.gradReadout);
            double grad_area_max = gradReadoutFlyback.getTotalAbsArea();
            double grad_area_sequence_max = 100 * (timeGradTopFlyback + grad_shape_rise_time);
            if (grad_area_max > grad_area_sequence_max) {
                double grad_time_flyback_min = parent.ceilToSubDecimal(grad_area_max / 100.0 - grad_shape_rise_time, 5);
                timeGradTopFlyback = grad_time_flyback_min;
                parent.getParam(UP.GRADIENT_FLYBACK_TIME).setValue(timeGradTopFlyback);
                parent.set(SP.Time_flyback, timeGradTopFlyback);
                gradReadoutFlyback.rePrepare();
            }
            gradReadoutFlyback.applyAmplitude();
        }
    }

    protected void setSeqParamTime() {
        timeGradTopFlyback = isFlyBackEnabled ? parent.getDouble(UP.GRADIENT_FLYBACK_TIME) : parent.minInstructionDelay;
        parent.set(SP.Time_flyback, timeGradTopFlyback);

        time_flyback_ramp = isFlyBackEnabled ? parent.getDouble(GRADIENT_RISE_TIME) : parent.minInstructionDelay;
        parent.set(SP.Time_grad_ramp_flyback, time_flyback_ramp);
    }
}
