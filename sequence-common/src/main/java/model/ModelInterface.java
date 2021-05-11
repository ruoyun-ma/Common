package model;

import common.RFPulse;
import kernel.SeqPrep;

public interface ModelInterface {
    void init(SeqPrep parent);
    void initFinal() throws Exception;
    void prep() throws Exception;
    void prepFinal();
    RFPulse getRfPulses();
    double getDuration();
    String getName();
    boolean isEnabled();
}
