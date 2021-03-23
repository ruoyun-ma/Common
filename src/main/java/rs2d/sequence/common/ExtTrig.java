package rs2d.sequence.common;

import rs2d.spinlab.sequence.SequenceTool;
import rs2d.spinlab.sequence.element.TimeElement;
import rs2d.spinlab.sequence.table.Table;
import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.Param;
import rs2d.spinlab.tools.param.TextParam;
import rs2d.spinlab.tools.table.Order;

import java.util.List;

import static java.util.Arrays.asList;

public class ExtTrig implements ModelInterface{
    private SeqPrep parent;
    protected static boolean isTriggerEnabled;

    public List<Double> triggerTime;
    public int nb_trigger;
    public double time_external_trigger_delay_max;
    public Table triggerdelay;
    public TextParam triggerChanel;

    protected enum UP implements GeneratorParamEnum {
        TRIGGER_EXTERNAL,
        TRIGGER_TIME,
        TRIGGER_CHANEL;

        @Override
        public Param build() {
            return null;
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum {
        Synchro_trigger,
        Ext_trig_source,
        Time_trigger_delay;
    }

    public ExtTrig(SeqPrep parent) {
        this.parent = parent;
    }

    @Override
    public void init() {
        triggerChanel = parent.getParam(UP.TRIGGER_CHANEL);
        triggerChanel.setSuggestedValues(asList(
                SequenceTool.ExtTrigSource.Ext1.name(),
                SequenceTool.ExtTrigSource.Ext2.name(),
                SequenceTool.ExtTrigSource.Ext1_AND_Ext2.name(),
                SequenceTool.ExtTrigSource.Ext1_XOR_Ext2.name()));
        triggerChanel.setRestrictedToSuggested(true);


        isTriggerEnabled = parent.getBoolean(UP.TRIGGER_EXTERNAL);
        triggerTime = parent.getListDouble(UP.TRIGGER_TIME);
        nb_trigger = isTriggerEnabled ? triggerTime.size() : 1;
        isTriggerEnabled = isTriggerEnabled && (nb_trigger > 0);
    }

    @Override
    public void initFinal() {}

    @Override
    public void prep() {
        parent.set(SP.Synchro_trigger, isTriggerEnabled ? TimeElement.Trigger.External : TimeElement.Trigger.Timer);
        parent.getSequenceParam(SP.Synchro_trigger).setLocked(true);
        time_external_trigger_delay_max = parent.minInstructionDelay;

        triggerdelay = parent.setSequenceTableValues(SP.Time_trigger_delay, Order.Four);

        if ((!isTriggerEnabled)) {
            triggerdelay.add(parent.minInstructionDelay);
        } else {
            for (int i = 0; i < nb_trigger; i++) {
                double time_external_trigger_delay = parent.roundToDecimal(triggerTime.get(i), 7);
                time_external_trigger_delay = Math.max(time_external_trigger_delay, parent.minInstructionDelay);
                triggerdelay.add(time_external_trigger_delay);
                time_external_trigger_delay_max = Math.max(time_external_trigger_delay_max, time_external_trigger_delay);
            }
        }
        parent.set(SP.Ext_trig_source, UP.TRIGGER_CHANEL);
    }

    @Override
    public void prepFinal() {
        if (parent.hasParam(CommonUP.SEQ_DESCRIPTION)) {
            String seqDescription = parent.getText(CommonUP.SEQ_DESCRIPTION);

            if (isTriggerEnabled && nb_trigger != 1) {
                seqDescription += "_TRIG=" + nb_trigger;
            } else if (isTriggerEnabled) {
                seqDescription += "_TRIG";
            }

            parent.getParam(CommonUP.SEQ_DESCRIPTION).setValue(seqDescription);
        }
    }

    @Override
    public boolean isEnabled() {
        return isTriggerEnabled;
    }

    @Override
    public RFPulse getRfPulses() {
        return null;
    }

    @Override
    public double getDuration() {
        return triggerdelay.getMaxValue();
    }

    @Override
    public String getName() {
        return "ExtTrig";
    }

}
