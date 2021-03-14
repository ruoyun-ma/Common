package rs2d.sequence.common;

import java.util.TreeMap;

public interface ModelInterface {
    void init() throws Exception;
    void initFinal() throws Exception;
    void prep() throws Exception;
    void prepFinal();
    RFPulse getRfPulses();
    boolean isEnabled();
}
