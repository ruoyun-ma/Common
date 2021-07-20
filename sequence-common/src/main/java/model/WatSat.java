package model;

import common.*;
import kernel.*;
import rs2d.commons.log.Log;
import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.NumberParam;
import rs2d.spinlab.tools.param.Param;
import rs2d.spinlab.tools.table.Order;

import java.util.Arrays;

import static common.CommonUP.*;
import static common.CommonSP.*;

/**
 * V1.0 - 2021.3 XG
 */


public class WatSat implements ModelInterface {
    public SeqPrep parent;
    boolean isWatSatEnabled;

    public RFPulse pulseTXWatSat;
    public Gradient gradFatsatRead;
    public Gradient gradFatsatPhase;
    public Gradient gradFatsatSlice;
    public int nb_watsat;

    private Order LoopOrder = Order.LoopB;
    protected boolean isAttAuto;
    protected double[] flipAngle;

    protected enum UP implements GeneratorParamEnum {
        WATSAT_ENABLED,
        WATSAT_TX_SHAPE,
        WATSAT_TX_LENGTH,
        WATSAT_BANDWIDTH,
        WATSAT_SP_FACTOR,
        WATSAT_FLIP_ANGLE,
        WATSAT_OFFSET_FREQ,
        WATSAT_GAMMA_B1,
        WATSAT_SATWAT,
        WATSAT_DELAY,
        ;

        @Override
        public Param build() {
            throw new UnsupportedOperationException("Unable to create a Param");
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum {
        Enable_ws,
        Tx_amp_ws,
        Tx_phase_ws,
        Tx_shape_ws,
        Tx_shape_phase_ws,
        Tx_att_offset_ws,
        Time_tx_ws,
        Time_grad_ws,
        Time_grad_ramp_ws,
        Time_delay_ws,
        Freq_offset_tx_ws,

        Grad_amp_ws_phase,
        Grad_amp_ws_read,
        Grad_amp_ws_slice,
        ;
    }

    public WatSat() {
    }

    @Override
    public void init(SeqPrep parent) {
        this.parent = parent;
        parent.setSuggestedValFromListString(parent.tx_shape, true, UP.WATSAT_TX_SHAPE);
    }

    @Override
    public void initPre() throws Exception {
        isWatSatEnabled = parent.getBoolean(UP.WATSAT_ENABLED);
        isAttAuto = parent.getBoolean(TX_AMP_ATT_AUTO);
    }

    @Override
    public void initFinal() throws Exception {
        parent.set(SP.Enable_ws, isWatSatEnabled);
        setSeqParamTime();

        if (isWatSatEnabled)
            nb_watsat = parent.getListInt(UP.WATSAT_FLIP_ANGLE).size();
        else
            nb_watsat = 1;

        flipAngle = new double[nb_watsat];
        if (false) {
            if (isWatSatEnabled)
                flipAngle = parent.getListDouble(UP.WATSAT_FLIP_ANGLE).stream().mapToDouble(d -> d).toArray();
            else
                flipAngle[0] = 0.0;
        } else {
            if (isWatSatEnabled) {
                flipAngle[0] = 360.0 * parent.getDouble(UP.WATSAT_TX_LENGTH) * parent.getDouble(UP.WATSAT_GAMMA_B1);
                parent.getParam(UP.WATSAT_FLIP_ANGLE).setValue(flipAngle[0]);
            } else
                flipAngle[0] = 0.0;
        }

        initPulseandGrad();
    }

    @Override
    public void prep() throws Exception {
        prepPulse();
        prepGrad();
        prepPulseComp();
    }

    @Override
    public void prepFinal() {
        Log.warning(getClass(), ("The order of watsat is LoopOrder = " + LoopOrder.name()));
        pulseTXWatSat.prepTxAmpMultiFA(parent.getListInt(TX_ROUTE), flipAngle, LoopOrder);
    }

    @Override
    public RFPulse getRfPulses() {
        return pulseTXWatSat;
    }

    @Override
    public String getName() {
        return "WatSat";
    }

    @Override
    public double getDuration() {
        double GradDuration = 2 * parent.getSequenceTable(SP.Time_grad_ramp_ws).get(0).doubleValue()
                + parent.getSequenceTable(SP.Time_grad_ws).get(0).doubleValue();
        double RFDuration = 2 * parent.getSequenceTable(SP.Time_grad_ramp_ws).get(0).doubleValue()
                + parent.getSequenceTable(SP.Time_tx_ws).get(0).doubleValue();
        double OverlapDuration = parent.getSequenceTable(SP.Time_grad_ramp_ws).get(0).doubleValue();
        double Delay = parent.getSequenceTable(SP.Time_delay_ws).get(0).doubleValue();
        double eachDuration = GradDuration + RFDuration - OverlapDuration + Delay;

        return eachDuration * nb_watsat;
    }

    @Override
    public boolean isEnabled() {
        return isWatSatEnabled;
    }

    protected void setSeqParamTime() {
        if (parent.hasParam(UP.WATSAT_DELAY))
            parent.set(SP.Time_delay_ws, UP.WATSAT_DELAY);
        else
            parent.set(SP.Time_delay_ws, parent.minInstructionDelay);

        double tx_bandwidth_90_ws = parent.getDouble(UP.WATSAT_BANDWIDTH);
        double tx_bandwidth_factor_90_ws = parent.getTx_bandwidth_factor(UP.WATSAT_TX_SHAPE, TX_BANDWIDTH_FACTOR, TX_BANDWIDTH_FACTOR_3D);
        double tx_length_90_ws = isWatSatEnabled ? tx_bandwidth_factor_90_ws / tx_bandwidth_90_ws : parent.minInstructionDelay;

        parent.getParam(UP.WATSAT_TX_LENGTH).setValue(tx_length_90_ws);
        parent.set(SP.Time_tx_ws, tx_length_90_ws);
    }

    protected void prepPulse() throws Exception {
        if (isAttAuto)
            clcPulse();
    }

    protected void prepGrad() {
        double pixmax;
        if (parent.hasParam(RESOLUTION_FREQUENCY) && parent.hasParam(RESOLUTION_PHASE) && parent.hasParam(RESOLUTION_SLICE)) {
            pixmax = Math.max(Math.max(parent.getDouble(RESOLUTION_FREQUENCY), parent.getDouble(RESOLUTION_PHASE)),
                    parent.getDouble(RESOLUTION_SLICE));
        } else {
            pixmax = Math.max(Math.max(parent.getDouble(FIELD_OF_VIEW), parent.getDouble(FIELD_OF_VIEW_PHASE)),
                    parent.getDouble(FIELD_OF_VIEW_3D));
        }

        if (isWatSatEnabled) {
            if (parent.hasParam(UP.WATSAT_SP_FACTOR)) {
                gradFatsatRead.setSpoiler(parent.getDouble(UP.WATSAT_SP_FACTOR), pixmax, 0, 0);
                gradFatsatPhase.setSpoiler(parent.getDouble(UP.WATSAT_SP_FACTOR), 0, pixmax, 0);
                gradFatsatSlice.setSpoiler(parent.getDouble(UP.WATSAT_SP_FACTOR), 0, 0, pixmax);
            } else {
                gradFatsatRead.setSpoiler(3.0, pixmax, 0, 0);
                gradFatsatPhase.setSpoiler(3.0, 0, pixmax, 0);
                gradFatsatSlice.setSpoiler(3.0, 0, 0, pixmax);
            }
        } else {
            gradFatsatRead.setSpoiler(0.0, pixmax, 0, 0);
            gradFatsatPhase.setSpoiler(0.0, 0, pixmax, 0);
            gradFatsatSlice.setSpoiler(0.0, 0, 0, pixmax);
        }

        gradFatsatRead.applyAmplitude(isWatSatEnabled ? LoopOrder : Order.FourLoop);
        gradFatsatPhase.applyAmplitude(isWatSatEnabled ? LoopOrder : Order.FourLoop);
        gradFatsatSlice.applyAmplitude(isWatSatEnabled ? LoopOrder : Order.FourLoop);
    }

    protected void prepPulseComp() {
        pulseTXWatSat.setFrequencyOffset(getFrequency() - parent.observeFrequency);
    }

    protected void initPulseandGrad() throws Exception {
        pulseTXWatSat = RFPulse.createRFPulse(parent.getSequence(), Tx_att, SP.Tx_amp_ws, SP.Tx_phase_ws,
                SP.Time_tx_ws, SP.Tx_shape_ws, SP.Tx_shape_phase_ws, SP.Freq_offset_tx_ws);
        if (parent.getSequence().getPublicTable(SP.Tx_att_offset_ws.name()) != null) {
            pulseTXWatSat.createAttOffset(parent.getSequence(), SP.Tx_att_offset_ws);
        }
        pulseTXWatSat.setShape(parent.getText(UP.WATSAT_TX_SHAPE), parent.nb_shape_points, "Hamming");

        gradFatsatRead = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_ws_read, SP.Time_grad_ws,
                Grad_shape_rise_up, Grad_shape_rise_down, SP.Time_grad_ramp_ws, parent.nucleus);
        gradFatsatPhase = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_ws_phase, SP.Time_grad_ws,
                Grad_shape_rise_up, Grad_shape_rise_down, SP.Time_grad_ramp_ws, parent.nucleus);
        gradFatsatSlice = Gradient.createGradient(parent.getSequence(), SP.Grad_amp_ws_slice, SP.Time_grad_ws,
                Grad_shape_rise_up, Grad_shape_rise_down, SP.Time_grad_ramp_ws, parent.nucleus);
    }

    protected void clcPulse() {
        if (!pulseTXWatSat.checkPower(getFlipAngle(), getFrequency(), parent.nucleus)) {
            parent.notifyOutOfRangeParam(UP.WATSAT_TX_LENGTH, pulseTXWatSat.getPulseDuration(),
                    ((NumberParam) parent.getParam(UP.WATSAT_TX_LENGTH)).getMaxValue(), "Pulse length too short to reach RF power with this pulse shape");
            parent.set(SP.Time_tx_ws, pulseTXWatSat.getPulseDuration());
            parent.getParam(UP.WATSAT_TX_LENGTH).setValue(pulseTXWatSat.getPulseDuration());
        }
    }

    protected double getFlipAngle() {
        if (flipAngle.length > 1)
            return Arrays.stream(flipAngle).max().getAsDouble();
        else
            return flipAngle[0];
    }

    protected double getFrequency() {
        if (parent.hasParam(UP.WATSAT_SATWAT)) {
            switch (parent.getText(UP.WATSAT_SATWAT)) {
                case "WATER":
                    return parent.observeFrequency;

                case "LIPID":
                    return parent.observeFrequency + parent.getDouble(FatSat.UP.FATSAT_OFFSET_FREQ);

                case "FREE":
                    return parent.observeFrequency + parent.getDouble(UP.WATSAT_OFFSET_FREQ);

                default:
                    return parent.observeFrequency + parent.getDouble(FatSat.UP.FATSAT_OFFSET_FREQ) / 2;
            }
        } else
            return parent.observeFrequency + parent.getDouble(FatSat.UP.FATSAT_OFFSET_FREQ) / 2;
    }


    public void setOrder(Order order) {
        // we provide an API for user input order
        LoopOrder = order;
    }
}