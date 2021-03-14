package rs2d.sequence.common;

import rs2d.spinlab.sequenceGenerator.GeneratorParamEnum;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.Param;


public class FatSatWep extends FatSat {
    private SeqPrep parent;

    public RFPulse pulseTXFatSat;

    //TODO: I really thought we should unify the name one day
    protected enum UP implements GeneratorParamEnum {
        FAT_SATURATION_WEP_ENABLED,
        FATSAT_WEP_TX_LENGTH;

        @Override
        public Param build() {
            return null;
        }
    }

    protected enum SP implements GeneratorSequenceParamEnum{
        Time_tx_fatsat_wep,
        Time_before_fatsat_wep_pulse,
        Tx_phase_fatsat_wep,
        Freq_offset_tx_fatsat_wep;
    }

    public FatSatWep(SeqPrep parent) {
        super(parent);
    }

    @Override
    public void init() throws Exception {
        super.init();
        isFatSatWepEnabled = parent.getBoolean(UP.FAT_SATURATION_WEP_ENABLED);

        if (isFatSatWepEnabled) { //FatSatWepEnabled has a higher priority than FatSatEnabled
            isFatSatEnabled = false;
            if (parent.hasParam(FatSat.UP.FAT_SATURATION_ENABLED))
                parent.getParam(FatSat.UP.FAT_SATURATION_ENABLED).setValue(isFatSatEnabled);
        }
    }

    @Override
    public boolean isEnabled() {
        return isFatSatWepEnabled;
    }

    @Override
    protected void setSeqParamTime() {
        super.setSeqParamTime();

        parent.set(SP.Time_before_fatsat_wep_pulse, parent.blankingDelay);
        if (isFatSatWepEnabled){
            parent.set(SP.Time_tx_fatsat_wep, parent.getDouble(UP.FATSAT_WEP_TX_LENGTH));
        }else{
            parent.set(SP.Time_tx_fatsat_wep, parent.minInstructionDelay);
        }
    }

}
