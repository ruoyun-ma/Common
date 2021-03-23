package rs2d.sequence.common;

public interface ModelInterface {
    void init();
    void initFinal() throws Exception;
    void prep() throws Exception;
    void prepFinal();
    RFPulse getRfPulses();
    double getDuration();
    String getName();
    boolean isEnabled();
}
