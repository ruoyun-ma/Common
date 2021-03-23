package rs2d.sequence.common;

import rs2d.spinlab.sequence.table.Table;
import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.Param;
import rs2d.spinlab.tools.table.Order;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 *  V1.0 - 2021.3 XG
 *
 */

public class InvRec implements ModelInterface {
    private SeqPrep parent;
    protected static boolean isInvRecEnabled;
    private double time0_IR_90 = 0.0;

    protected RFPulse pulseTXIR;

    public List<Double> inversionRecoveryTime;
    public int nb_inversionRecovery;
    public double time_IR_delay_max;
    public Table time_TI_delay;

    protected enum UP implements GeneratorParamEnum {
        INVERSION_RECOVERY_ENABLED,
        INVERSION_TIME_MULTI,
        GRADIENT_ENABLE_CRUSHER_IR,
        GRADIENT_CRUSHER_IR_TOP_TIME,
        GRADIENT_AMP_CRUSHER_IR,
        INVERSION_RECOVERY_TAU, // dedicated input from the user
        ;

        @Override
        public Param build() {
            return null;
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum {
        Grad_enable_IR,
        Tx_enable_IR,
        Grad_enable_crush_IR,
        Grad_amp_crusher_IR,
        Time_grad_IR_crusher_top,
        Time_grad_IR_crusher_ramp,
        Time_tx_IR_length,
        Time_grad_IR_ramp,
        Time_TI_delay,
        Tx_freq_offset_IR,
        FreqOffset_tx_prep_IR,
        ;
    }

    public InvRec(SeqPrep parent) {
        this.parent = parent;
    }

    @Override
    public void init() {
        time_TI_delay = parent.getSequenceTable(SP.Time_TI_delay);

        isInvRecEnabled = parent.getBoolean(UP.INVERSION_RECOVERY_ENABLED);
        inversionRecoveryTime = parent.getListDouble(UP.INVERSION_TIME_MULTI);
        nb_inversionRecovery = isInvRecEnabled ? inversionRecoveryTime.size() : 1;
        isInvRecEnabled = isInvRecEnabled && (nb_inversionRecovery >= 1);
    }

    @Override
    public void initFinal() throws Exception {
        // Avoid multi IR time when Multi echo
        if (nb_inversionRecovery != 1 && (parent.nb_echo_4D != 1)) {
            double tmp = inversionRecoveryTime.get(0);
            inversionRecoveryTime.clear();
            inversionRecoveryTime.add(tmp);
            nb_inversionRecovery = 1;
        }

        parent.set(SP.Grad_enable_IR, isInvRecEnabled);
        parent.set(SP.Tx_enable_IR, isInvRecEnabled);
        parent.set(SP.Grad_enable_crush_IR, isInvRecEnabled && parent.getBoolean(UP.GRADIENT_ENABLE_CRUSHER_IR));


        setSeqParamTime();
        initPulseandGrad();
    }

    @Override
    public void prep() throws Exception {
        prepPulseComp();
    }

    @Override
    public void prepFinal() {
        // ------------------------------------------
        // calculate delays adapted to current TI
        // ------------------------------------------

        time_IR_delay_max = 0.00001;
        time_TI_delay.clear();

        //  double time_TI_delay;
        if (isInvRecEnabled) {
            //TI
            if (parent.hasParam(UP.INVERSION_RECOVERY_TAU) && parent.getDouble(UP.INVERSION_RECOVERY_TAU) > 0.0) {
                time0_IR_90 = parent.getDouble(UP.INVERSION_RECOVERY_TAU);
            } else {
                time0_IR_90 += parent.getSequenceTable(SP.Time_tx_IR_length).get(0).doubleValue() / 2
                        + parent.getSequenceTable(SP.Time_grad_IR_ramp).get(0).doubleValue()
                        + 2 * parent.getSequenceTable(SP.Time_grad_IR_crusher_ramp).get(0).doubleValue()
                        + parent.getSequenceTable(SP.Time_grad_IR_crusher_top).get(0).doubleValue()
                        + parent.minInstructionDelay;

                time0_IR_90 += parent.getDouble(CommonUP.GRADIENT_RISE_TIME)
                        + (parent.hasParam(CommonUP.TX_LENGTH_90) ? parent.getDouble(CommonUP.TX_LENGTH_90) : parent.getDouble(CommonUP.TX_LENGTH)) / 2;
            }
            if (parent.models.containsKey("SatBand")) {
                time0_IR_90 += parent.models.get("SatBand").getDuration() + parent.minInstructionDelay;
            }
            if (parent.models.containsKey("FatSat") && parent.models.containsKey("FatSatWep")) {
                time0_IR_90 += parent.models.get("FatSatWep").getDuration();
            } else if (parent.models.containsKey("FatSat")) {
                time0_IR_90 += parent.models.get("FatSat").getDuration();
            }

            boolean increaseTI = false;
            // ArrayList<Number> arrayListTI = new ArrayList<Number>();
            Collection<Double> arrayListTI_min = new ArrayList<Double>();
            for (int i = 0; i < nb_inversionRecovery; i++) {
                double IR_time = inversionRecoveryTime.get(i);
                // arrayListTI.add(IR_time);
                double delay0 = IR_time - time0_IR_90;
                System.out.println("delay0 " + delay0);

                if ((delay0 < parent.minInstructionDelay)) {
                    double ti_min = parent.ceilToSubDecimal((time0_IR_90 + parent.minInstructionDelay), 4);
                    System.out.println("IR_time-" + IR_time + " -- >" + ti_min + " car " + delay0);
                    IR_time = ti_min;
                    increaseTI = true;
                }
                arrayListTI_min.add(IR_time);
                time_IR_delay_max = Math.max(time_IR_delay_max, IR_time);
                delay0 = IR_time - time0_IR_90;
                time_TI_delay.add(delay0);
            }


            if (increaseTI) {
                System.out.println(" --- -- -- - - - increaseTI-------------------- ");

                inversionRecoveryTime.clear();
                inversionRecoveryTime.addAll(arrayListTI_min);
                //this.notifyOutOfRangeParam(rs2d.sequence.spinecho.SPIN_ECHO_devParams.INVERSION_TIME_MULTI.name(), arrayListTI, arrayListTI_min, ((NumberParam) getParam(rs2d.sequence.spinecho.SPIN_ECHO_devParams.INVERSION_TIME_MULTI")).getMaxValue(), "TI too short");
                parent.getParam(UP.INVERSION_TIME_MULTI).setValue(inversionRecoveryTime);
            }
            time_TI_delay.setOrder(Order.Four);
            time_TI_delay.setLocked(true);
        } else {
            time_IR_delay_max = parent.minInstructionDelay;
            time_TI_delay.add(parent.minInstructionDelay);
        }


        // ------------------------------------------
        // SEQ_DESCRIPTION
        // ------------------------------------------
        if (parent.hasParam(CommonUP.SEQ_DESCRIPTION)) {
            String seqDescription = parent.getText(CommonUP.SEQ_DESCRIPTION);

            if (isInvRecEnabled && nb_inversionRecovery != 1) {
                seqDescription += "_IR-" + nb_inversionRecovery;
            } else if (isInvRecEnabled) {
                seqDescription += "_IR=" + inversionRecoveryTime.get(0) + "s";
            }

            parent.getParam(CommonUP.SEQ_DESCRIPTION).setValue(seqDescription);
        }
    }

    @Override
    public RFPulse getRfPulses() {
        return null;
    }

    @Override
    public double getDuration() {
        return time0_IR_90 + time_IR_delay_max;
    }

    @Override
    public String getName() {
        return "InvRec";
    }

    @Override
    public boolean isEnabled() {
        return isInvRecEnabled;
    }

    protected void initPulseandGrad() {
        pulseTXIR = RFPulse.createRFPulse(parent.getSequence(), CommonSP.Tx_att, CommonSP.Tx_amp_180, CommonSP.Tx_phase_180,
                SP.Time_tx_IR_length, CommonSP.Tx_shape_180, CommonSP.Tx_shape_phase_180, SP.Tx_freq_offset_IR);
    }

    protected void prepPulseComp() {
        if (parent.isMultiplanar && parent.nb_planar_excitation > 1 && parent.isEnableSlice) {
            if (isInvRecEnabled) {
                pulseTXIR.prepareOffsetFreqMultiSlice(parent.gradSlice180, parent.nb_planar_excitation, parent.spacingBetweenSlice, parent.off_center_distance_3D);
                //TODO:XG: we need to check if it works for VFL cases
                pulseTXIR.reoderOffsetFreq(parent.plugin, parent.acqMatrixDimension1D * parent.echoTrainLength, parent.nb_slices_acquired_in_single_scan);
            }
            pulseTXIR.setFrequencyOffset(isInvRecEnabled ? (parent.nb_slices_acquired_in_single_scan != 1 ? Order.ThreeLoop : Order.Three) : Order.FourLoop);
        } else {
            if (isInvRecEnabled) {
                pulseTXIR.prepareOffsetFreqMultiSlice(parent.gradSlice180, 1, 0, parent.off_center_distance_3D);
            }
            pulseTXIR.setFrequencyOffset(isInvRecEnabled ? Order.Three : Order.FourLoop);
        }

        RFPulse pulseTXIRPrep = RFPulse.createRFPulse(parent.getSequence(), CommonSP.Time_grad_ramp, SP.FreqOffset_tx_prep_IR);
        pulseTXIRPrep.setCompensationFrequencyOffset(pulseTXIR, 0.5);
    }

    protected void setSeqParamTime() {
        if (isInvRecEnabled) {
            boolean is_IR_crusher = parent.getBoolean(UP.GRADIENT_ENABLE_CRUSHER_IR);
            // prepare Crusher gradient
            if (is_IR_crusher) {
                double time_grad_IR_tmp_top = parent.getDouble(UP.GRADIENT_CRUSHER_IR_TOP_TIME);

                // IR Crusher:
                double grad_amp_crusher_IR = parent.getDouble(UP.GRADIENT_AMP_CRUSHER_IR);
                parent.set(SP.Grad_amp_crusher_IR, grad_amp_crusher_IR);
                parent.set(SP.Time_grad_IR_crusher_top, time_grad_IR_tmp_top);
                parent.set(SP.Time_grad_IR_crusher_ramp, parent.getDouble(CommonUP.GRADIENT_RISE_TIME));
            } else {
                parent.set(SP.Time_grad_IR_crusher_top, parent.minInstructionDelay);
                parent.set(SP.Time_grad_IR_crusher_ramp, parent.minInstructionDelay);
            }
            parent.set(SP.Time_tx_IR_length, parent.getDouble(CommonUP.TX_LENGTH_180));
            parent.set(SP.Time_grad_IR_ramp, parent.getDouble(CommonUP.GRADIENT_RISE_TIME));
        } else {
            parent.set(SP.Time_tx_IR_length, parent.minInstructionDelay);
            parent.set(SP.Time_grad_IR_ramp, parent.minInstructionDelay);
            parent.set(SP.Time_grad_IR_crusher_top, parent.minInstructionDelay);
            parent.set(SP.Time_grad_IR_crusher_ramp, parent.minInstructionDelay);
        }
    }
}
