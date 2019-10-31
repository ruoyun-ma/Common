package rs2d.sequence.common;

import rs2d.spinlab.data.transformPlugin.TransformPlugin;
import rs2d.spinlab.instrument.Instrument;
import rs2d.spinlab.instrument.InstrumentTxChannel;
import rs2d.spinlab.instrument.probe.Probe;
import rs2d.spinlab.instrument.probe.ProbeChannelPower;
import rs2d.spinlab.instrument.util.GradientMath;
import rs2d.spinlab.instrument.util.TxMath;
import rs2d.spinlab.plugins.loaders.LoaderFactory;
import rs2d.spinlab.plugins.loaders.PluginLoaderInterface;
import rs2d.spinlab.sequence.Sequence;
import rs2d.spinlab.sequence.table.Shape;
import rs2d.spinlab.sequence.table.Table;
import rs2d.spinlab.sequence.table.Utility;
import rs2d.spinlab.sequence.table.generator.TableGeneratorInterface;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.NumberParam;
import rs2d.spinlab.tools.param.Param;
import rs2d.spinlab.tools.table.Order;
import rs2d.spinlab.tools.utility.Nucleus;

import java.util.List;

/**
 * Class RFPulse
 * V2.6 constructor with generatorSequenceParam .name() V2019.06
 * V2.5- getNearestSW Sup Inf for Cam4
 * V2.3- 2019-06-06 JR
 * V2.2- 2018-12-19 JR
 * V2.1- 2017-10-24 JR
 */
public class RFPulse {
    private Table amplitudeTable = null;
    private Param attParam = null;
    private Table phase = null;
    private Table FrequencyOffsetTable = null;

    private Table timeTable = null;
    private Shape shape = null;
    private Shape shapePhase = null;

    int numberOfFreqOffset = -1;
    Order FrequencyOffsetOrder = Order.FourLoopB;

    boolean isSlr = false;
    int slrIndex = -1;

    double[] slrPowerFactors90 = {2.331, 2.331 * 0.72}; //slr power factors compared to not slr pulses
    double[] slrPowerFactors180 = {3.879, 3.879 * 0.36};

    private double[] txFrequencyOffsetTable = null;
    private int txAtt = -1;

    private double tx_amp90 = Double.NaN;
    private double tx_amp180 = Double.NaN;
    private double tx_amp = Double.NaN;
    private double flipAngle = Double.NaN;

    private double pulseDuration;
    private double observeFrequency = Double.NaN;
    private Nucleus nucleus;
    private double powerPulse = Double.NaN;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  Constructor
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public RFPulse(Param attPar, Table amplitudeTab, Table phaseTab,
                   Table timeTab, Shape shape, Shape shapePhase, Table offsetFreqTab) {
        amplitudeTable = amplitudeTab;
        attParam = attPar;
        phase = phaseTab;

        timeTable = timeTab;
        this.shape = shape;
        this.shapePhase = shapePhase;

        pulseDuration = timeTable.get(0).doubleValue();
        FrequencyOffsetTable = offsetFreqTab;
    }


    public RFPulse(Table timeTab, Table offsetFreqTab) {
        timeTable = timeTab;

        pulseDuration = timeTab.get(0).doubleValue();
        FrequencyOffsetTable = offsetFreqTab;
    }

    public RFPulse(Table timeTab, Table offsetFreqTab, Table phaseTab) {
        timeTable = timeTab;
        phase = phaseTab;
        pulseDuration = timeTab.get(0).doubleValue();
        FrequencyOffsetTable = offsetFreqTab;
    }


    public static RFPulse createRFPulse(Sequence sequence, GeneratorSequenceParamEnum txAttParam, GeneratorSequenceParamEnum amplitudeTab, GeneratorSequenceParamEnum txPhaseTab,
                                        GeneratorSequenceParamEnum timeTab, GeneratorSequenceParamEnum shape, GeneratorSequenceParamEnum shapePhaseShape, GeneratorSequenceParamEnum offsetTab) {
        return new RFPulse(sequence.getPublicParam(txAttParam.name()), sequence.getTable(amplitudeTab.name()), sequence.getTable(txPhaseTab.name()),
                sequence.getPublicTable(timeTab.name()), (Shape) sequence.getTable(shape.name()), (Shape) sequence.getTable(shapePhaseShape.name()), sequence.getTable(offsetTab.name()));
    }

    public static RFPulse createRFPulse(Sequence sequence, GeneratorSequenceParamEnum timeTab, GeneratorSequenceParamEnum offsetTab) {
        return new RFPulse(sequence.getPublicTable(timeTab.name()), sequence.getTable(offsetTab.name()));
    }

    public static RFPulse createRFPulse(Sequence sequence, GeneratorSequenceParamEnum timeTab, GeneratorSequenceParamEnum offsetTab, GeneratorSequenceParamEnum txPhaseTab) {
        return new RFPulse(sequence.getPublicTable(timeTab.name()), sequence.getTable(offsetTab.name()), sequence.getTable(txPhaseTab.name()));
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  general  methods: get set
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public double getPulseDuration() {
        return pulseDuration;
    }

    public int getAtt() {
        return txAtt;
    }

    public double getAmp() {
        return tx_amp;
    }

    public double getAmp90() {
        return tx_amp90;
    }

    public double getAmp180() {
        return tx_amp180;
    }

    public Order getFrequencyOffsetOrder() {
        return FrequencyOffsetOrder;
    }

    public int getNumberOfFreqOffset() {
        return numberOfFreqOffset;
    }

    public boolean isSlr() {
        return isSlr;
    }

    public double getSlrPowerFactors90() {
        return isSlr ? slrPowerFactors90[slrIndex] : 1.0;
    }

    public double getSlrPowerFactors180() {
        return isSlr ? slrPowerFactors180[slrIndex] : 1.0;
    }

    public void setAtt(int att) {
        txAtt = att;
        attParam.setValue(txAtt);

    }

    public void setAtt(NumberParam att) {
        txAtt = att.getValue().intValue();
        attParam.setValue(txAtt);

    }

    public void setAmp(double amp) {
        tx_amp = amp;
        setSequenceTableSingleValue(amplitudeTable, tx_amp);
        flipAngle = tx_amp * 90 / tx_amp90;
    }

    public void setAmp(Order order, double... amps) {
        setSequenceTableValues(amplitudeTable, order, amps);
        if (amps.length > 0) {
            tx_amp = amps[0];
        }
    }

    public void setAmp(NumberParam ampParam) {
        tx_amp = ampParam.getValue().doubleValue();
        setSequenceTableSingleValue(amplitudeTable, tx_amp);
        flipAngle = tx_amp * 90 / tx_amp90;
    }

    public void setAmp90(double amp90) {
        tx_amp90 = amp90;
    }

    public void setAmp180(double amp180) {
        tx_amp180 = amp180;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  Calculation Amp Att...
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     * Automatic Calibration for a given Flip angle with ATT set to get 180° -> 80%amp.
     * set amplitudeTable and attParam
     *
     * @param flipAngle         : flip angle of the pulse
     * @param observe_frequency :set pulse property
     * @param txRoute
     * @param nucleus           :set pulse propert
     * @return test_change_time : false if the pulse duration has changed
     */
    public boolean setAutoCalibFor180(double flipAngle, double observe_frequency, List<Integer> txRoute, Nucleus nucleus) {
        InstrumentTxChannel txCh = Instrument.instance().getTxChannels().get(txRoute.get(0));
        this.flipAngle = flipAngle;
        observeFrequency = observe_frequency;
        this.nucleus = nucleus;
        boolean test_change_time = true;
        if (txAtt == -1) {
            test_change_time = prepTxAttFor180(80, txCh);
        }
        // Calculate amp with the new and real Att
        tx_amp90 = calculateTxAmp90(txCh);
        tx_amp180 = attParamTxAmp180(txCh);
        tx_amp = tx_amp90 * flipAngle / 90;
        setSequenceTableSingleValue(amplitudeTable, tx_amp);
        attParam.setValue(txAtt);
        // set calculated parameters to display values & sequence
        return test_change_time;
    }

    /**
     * Automatic Calibration: ATT set to get 180°-> tx_amp%.
     *
     * @param tx_amp : of the 180° pulse
     * @param txCh
     * @return test_change_time : false if the pulse duration has changed
     */
    private boolean prepTxAttFor180(double tx_amp, InstrumentTxChannel txCh) {
        boolean test_change_time = true;
        Probe probe = Instrument.instance().getTransmitProbe();
        ProbeChannelPower pulse = TxMath.getProbePower(probe, null, nucleus.name());
        double instrument_length_180 = pulse.getHardPulse180().x;
        double power_factor = Utility.powerFillingFactor(shape) / getSlrPowerFactors180();       // get RF pulse power factor from instrument to calculate RF pulse amplitude
        // in future version V2019.07: double power_factor = Utility.complexPowerFillingFactor(shape) / getSlrPowerFactors180();       // get RF pulse power factor from instrument to calculate RF pulse amplitude
        double instrument_power_180 = pulse.getHardPulse180().y / power_factor;
        double power_180 = instrument_power_180 * Math.pow(instrument_length_180 / pulseDuration, 2);

        if (power_180 > pulse.getMaxRfPowerPulsed()) {  // TX LENGTH 90 MIN
            pulseDuration = ceilToSubDecimal(instrument_length_180 / Math.sqrt(pulse.getMaxRfPowerPulsed() / instrument_power_180), 6);
            setSequenceTableSingleValue(timeTable, pulseDuration);
            power_180 = instrument_power_180 * Math.pow(instrument_length_180 / pulseDuration, 2);
            test_change_time = false;
        }
        // Calculate Att to get a 180° RF pulse around 80% amp
        if (txAtt == -1)
            txAtt = 1 + (int) TxMath.getTxAttFor(power_180, txCh, tx_amp, observeFrequency);
        return test_change_time;
    }

    /**
     * calculate 90° Pulse Amplitude for txAtt
     *
     * @param txCh
     * @return tx_amp
     */
    private double calculateTxAmp90(InstrumentTxChannel txCh) {
        if (txAtt == -1) {
            txAtt = (int) attParam.getValue();
        }
        double tx_amp;
        Probe probe = Instrument.instance().getTransmitProbe();
        ProbeChannelPower pulse = TxMath.getProbePower(probe, null, nucleus.name());
        double power_factor = Utility.powerFillingFactor(shape) / getSlrPowerFactors90();       // get RF pulse power factor from instrument to calculate RF pulse amplitude
        // in future version V2019.07: double power_factor = Utility.complexPowerFillingFactor(shape) / getSlrPowerFactors180();       // get RF pulse power factor from instrument to calculate RF pulse amplitude
        double instrument_power_90 = pulse.getHardPulse90().y / power_factor;
        double instrument_length = pulse.getHardPulse90().x;
        double power = instrument_power_90 * Math.pow(instrument_length / pulseDuration, 2);
        tx_amp = TxMath.getTxAmpFor(power, txCh, txAtt, observeFrequency);
        return tx_amp;
    }

    /**
     * calculate 180° Pulse Amplitude for txAtt
     *
     * @param txCh
     * @return tx_amp
     */
    private double attParamTxAmp180(InstrumentTxChannel txCh) {
        double tx_amp;
        Probe probe = Instrument.instance().getTransmitProbe();
        ProbeChannelPower pulse = TxMath.getProbePower(probe, null, nucleus.name());
        double power_factor = Utility.powerFillingFactor(shape) / getSlrPowerFactors180();       // get RF pulse power factor from instrument to calculate RF pulse amplitude
        // in future version V2019.07: double power_factor = Utility.complexPowerFillingFactor(shape) ;       // get RF pulse power factor from instrument to calculate RF pulse amplitude
        double instrument_power_180 = pulse.getHardPulse180().y / power_factor;
        double instrument_length = pulse.getHardPulse180().x;
        double power = instrument_power_180 * Math.pow(instrument_length / pulseDuration, 2);
        tx_amp = TxMath.getTxAmpFor(power, txCh, txAtt, observeFrequency);
        return tx_amp;
    }

    /**
     * calculate 180° Pulse Amplitude for txAtt
     *
     * @param flipAngle         set pulse property
     * @param observe_frequency set pulse property
     * @param nucleus           set pulse property
     * @return test_change_time = false if the pulseDuration was increase because of exceeding power max
     */
    public boolean checkPower(double flipAngle, double observe_frequency, Nucleus nucleus) {
        this.flipAngle = flipAngle;
        observeFrequency = observe_frequency;
        this.nucleus = nucleus;
        boolean test_change_time = true;
        if (txAtt == -1) {
            test_change_time = calculatePower();
        }
        return test_change_time;
    }

    public double getPower() {
        if (Double.isNaN(powerPulse))
            calculatePower();
        return powerPulse;
    }

    /**
     * calculate powerPulse from flipAngle and pulseDuration
     *
     * @return test_change_time = false if the pulseDuration was increase because of exceeding power max
     */
    private boolean calculatePower() {
        if (Double.isNaN(flipAngle) || Double.isNaN(pulseDuration) || Double.isNaN(pulseDuration)) {
            System.out.println("Power Calculation cannot be done: Flip Angle or Time not set! ");
        }
        boolean test_change_time = true;
        Probe probe = Instrument.instance().getTransmitProbe();
        ProbeChannelPower pulse = TxMath.getProbePower(probe, null, nucleus.name());
        double power_factor = Utility.powerFillingFactor(shape) / (flipAngle < 135 ? getSlrPowerFactors90() : getSlrPowerFactors180());       // get RF pulse power factor from instrument to calculate RF pulse amplitude
        // in future version V2019.07: double power_factor = Utility.complexPowerFillingFactor(shape) ;       // get RF pulse power factor from instrument to calculate RF pulse amplitude
        double instrument_length = flipAngle < 135 ? pulse.getHardPulse90().x : pulse.getHardPulse180().x;
        double instrument_power = (flipAngle < 135 ? pulse.getHardPulse90().y : pulse.getHardPulse180().y) / power_factor;
        powerPulse = instrument_power * Math.pow(instrument_length / pulseDuration, 2) * Math.pow(flipAngle / (flipAngle < 135 ? 90 : 180), 2);
        if (powerPulse > pulse.getMaxRfPowerPulsed()) {  // TX LENGTH 90 MIN
            pulseDuration = ceilToSubDecimal(instrument_length / Math.sqrt(pulse.getMaxRfPowerPulsed() / (instrument_power * Math.pow(flipAngle / (flipAngle < 135 ? 90 : 180), 2))), 6);
            setSequenceTableSingleValue(timeTable, pulseDuration);
            powerPulse = instrument_power * Math.pow(instrument_length / pulseDuration, 2) * Math.pow(flipAngle / (flipAngle < 135 ? 90 : 180), 2);
            test_change_time = false;
        }
        // Calculate Att to get a 180° RF pulse around 80% amp
        return test_change_time;
    }

    /**
     * calculate and set txAtt in order to get the powerPulse at pulse amplitude = amp
     *
     * @param amp     amplitude for powerPulse
     * @param txRoute
     * @return txAtt
     */
    public int prepAtt(double amp, List<Integer> txRoute) {
        InstrumentTxChannel txCh = Instrument.instance().getTxChannels().get(txRoute.get(0));

        double tx_amp_180_desired = amp;     // set 180° RF puse arround 80% of Chameleon output
        txAtt = 1 + (int) TxMath.getTxAttFor(powerPulse, txCh, tx_amp_180_desired, observeFrequency);
        attParam.setValue(txAtt);
        return txAtt;
    }

    /**
     * prepare and set txAmp according to txAtt , flipAngle
     *
     * @param txRoute
     * @return tx_amp
     */
    public double prepTxAmp(List<Integer> txRoute) {
        if (txAtt == -1) {
            txAtt = ((NumberParam) attParam).getValue().intValue();
        }
        InstrumentTxChannel txCh = Instrument.instance().getTxChannels().get(txRoute.get(0));
        tx_amp90 = calculateTxAmp90(txCh);
        tx_amp180 = attParamTxAmp180(txCh);
        tx_amp = (flipAngle < 135 ? tx_amp90 : tx_amp180) * flipAngle / (flipAngle < 135 ? 90 : 180);
        setSequenceTableSingleValue(amplitudeTable, tx_amp);
        return tx_amp;
    }

    /**
     * prepare and set txAmp according to txAtt , flipAngle
     *
     * @param txRoute
     * @return tx_amp
     */
    public double prepTxAmpMultiFA(List<Integer> txRoute, double[] FA_list, Order order) {
        if (txAtt == -1) {
            txAtt = ((NumberParam) attParam).getValue().intValue();
        }
        InstrumentTxChannel txCh = Instrument.instance().getTxChannels().get(txRoute.get(0));
        tx_amp90 = calculateTxAmp90(txCh);
        tx_amp180 = attParamTxAmp180(txCh);
        setSequenceTableValues(amplitudeTable, order);
        for (int i = 0; i < FA_list.length; i++) {
            tx_amp = (FA_list[i] < 135 ? tx_amp90 : tx_amp180) * FA_list[i] / (FA_list[i] < 135 ? 90 : 180);
            amplitudeTable.add(tx_amp);
        }
        return tx_amp;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  Shape
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     * prepare and set txAmp according to txAtt , flipAngle
     *
     * @param pulseName     :GAUSSIAN , SINC3 , SINC5 , HARD
     * @param numberOfPoint The number of point of the generated shape
     * @param window
     */
    public void setShape(String pulseName, int numberOfPoint, String window) throws Exception {
        shapePhase.clear();
        if ("GAUSSIAN".equalsIgnoreCase(pulseName)) {
            setShapeTableValuesFromGaussGen(shape, numberOfPoint, 0.250, 100, false);
        } else if ("SINC3".equalsIgnoreCase(pulseName)) {
            setShapeTableValuesFromSincGen(shape, numberOfPoint, 2, 100, true, window);
            for (int i = 0; i < (numberOfPoint); i++) {
                shapePhase.add((i < ((int) (numberOfPoint / 4.0)))
                        || (i >= ((int) (numberOfPoint * 3.0 / 4.0))) ? 180 : 0);
            }
        } else if ("SINC5".equalsIgnoreCase(pulseName)) {
            setShapeTableValuesFromSincGen(shape, numberOfPoint, 3, 100, true, window);
            for (int i = 0; i < (numberOfPoint); i++) {
                shapePhase.add(((i > ((int) (numberOfPoint * 1.0 / 6.0)) && (i <= ((int) (numberOfPoint * 2.0 / 6.0))))
                        || (i > ((int) (numberOfPoint * 4.0 / 6.0)) && (i <= ((int) (numberOfPoint * 5.0 / 6.0)))))
                        ? 180 : 0);
            }

        } else if ("HARD".equalsIgnoreCase(pulseName)) {
            shape.clear();    // set to HARD pulse
            shapePhase.clear();    // set to HARD pulse
            shape.setFirst(100);
            shapePhase.setFirst(0);

        } else if ("RAMP".equalsIgnoreCase(pulseName)) {

            setTableValuesFromSincGenRamp(shape, numberOfPoint, 3, 100, true, window, 0.2);
            setTableValuesFromSincGenRampPhase(shapePhase, numberOfPoint, 3, 100, true, window, 0.2);
        } else if ("SLR_8_5152".equalsIgnoreCase(pulseName)) {
            shape.clear();
            shapePhase.clear();
            String bw = "8.5152";
            String type = "90 degree";
            setTableValuesFromSLRGen(shape, numberOfPoint, 0, 0, type, true, bw);
            setTableValuesFromSLRPhaseGen(shapePhase, numberOfPoint, 0, 0, type, false, bw);
            isSlr = true;
            slrIndex = 0;
        } else if ("SLR_4_2576".equalsIgnoreCase(pulseName)) {
            shape.clear();
            shapePhase.clear();
            String bw = "4.2576";
            String type = "90 degree";
            //type="Refocusing (spin-echo)";
            setTableValuesFromSLRGen(shape, numberOfPoint, 0, 0, type, true, bw);
            setTableValuesFromSLRPhaseGen(shapePhase, numberOfPoint, 0, 0, type, false, bw);
            isSlr = true;
            slrIndex = 1;
        }
    }

    private void setTableValuesFromSLRGen(Table table, int nbpoint, int bw, double amp, String type, boolean abs, String bwstring) throws Exception {
        TableGeneratorInterface gen = null;
        gen = loadTableGenerator("SLR");
        gen.getParams().get(0).setValue(nbpoint);
        gen.getParams().get(1).setValue(bwstring);
        gen.getParams().get(2).setValue(amp);
        gen.getParams().get(3).setValue(type);//type
        gen.getParams().get(4).setValue(abs);//abs

        table.setGenerator(gen);
        if (gen == null) {
            table.clear();
            table.setFirst(100);
        } else {
            gen.generate();
        }
    }

    private void setTableValuesFromSLRPhaseGen(Table table, int nbpoint, int bw, double amp, String type, boolean abs, String bwstring) throws Exception {
        TableGeneratorInterface gen = null;
        gen = loadTableGenerator("SLRPhase");
        gen.getParams().get(0).setValue(nbpoint);
        gen.getParams().get(1).setValue(bwstring);
        gen.getParams().get(2).setValue(amp);
        gen.getParams().get(3).setValue(type);//type
        gen.getParams().get(4).setValue(abs);//abs

        table.setGenerator(gen);
        if (gen == null) {
            table.clear();
            table.setFirst(100);
        } else {
            gen.generate();
        }
    }

    private void setTableValuesFromSincGenRamp(Table table, int nbpoint, int nblobe, double amp, Boolean abs, String window, double slope) throws Exception {
        TableGeneratorInterface gen = null;
        gen = loadTableGenerator("SincRamp");
        gen.getParams().get(0).setValue(nbpoint);
        gen.getParams().get(1).setValue(nblobe);
        gen.getParams().get(2).setValue(amp);
        gen.getParams().get(3).setValue(abs);//abs
        gen.getParams().get(4).setValue(window);//abs
        gen.getParams().get(5).setValue(slope);

        table.setGenerator(gen);
        if (gen == null) {
            table.clear();
            table.setFirst(100);
        } else {
            gen.generate();
        }
    }


    private void setTableValuesFromSincGenRampPhase(Table table, int nbpoint, int nblobe, double amp, Boolean abs, String window, double slope) throws Exception {
        TableGeneratorInterface gen = null;
        gen = loadTableGenerator("SincRampPhase");
        gen.getParams().get(0).setValue(nbpoint);
        gen.getParams().get(1).setValue(nblobe);
        gen.getParams().get(2).setValue(amp);
        gen.getParams().get(3).setValue(abs);//abs
        gen.getParams().get(4).setValue(window);//abs
        gen.getParams().get(5).setValue(slope);

        table.setGenerator(gen);
        if (gen == null) {
            table.clear();
            table.setFirst(100);
        } else {
            gen.generate();
        }
    }

    private TableGeneratorInterface loadTableGenerator(String generatorName) throws Exception {
        TableGeneratorInterface gen = null;
        PluginLoaderInterface loader = LoaderFactory.getTableGeneratorPluginLoader();
        if (loader.containsPlugin(generatorName)) {
            gen = (TableGeneratorInterface) loader.getPluginByName(generatorName);
        }
        return gen;
    }

    /**
     * Generate a table of element with a Sinus Cardinal generator
     *
     * @param table   The table to be set
     * @param nbpoint The number of point of the generated sinus cardinal
     * @param nblobe  The number of lob of the generated sinus cardinal
     * @param amp     The amplitude of the generated sinus cardinal (in %)
     * @param abs     true if you want the absolute values and false otherwise
     */
    private void setShapeTableValuesFromSincGen(Table table, int nbpoint, int nblobe, double amp, Boolean abs, String window) throws Exception {
        String name = "Sinus Cardinal with Apodisation";
        TableGeneratorInterface gen;
        gen = LoaderFactory.getTableGeneratorPluginLoader().getPluginByName(name);
        if (gen == null) {
            throw new IllegalStateException("Table generator not found: " + name);
        }
        table.setGenerator(gen);
        gen.getParams().get(0).setValue(nbpoint);
        gen.getParams().get(1).setValue(nblobe);
        gen.getParams().get(2).setValue(amp);
        gen.getParams().get(3).setValue(abs);//abs
        gen.getParams().get(4).setValue(window);//abs
        table.setGenerator(gen);
        gen.generate();
    }

    /**
     * Generate a table of element with a Gaussian generator
     *
     * @param table   The table to be set
     * @param nbpoint The number of point of the generated gaussian
     * @param width   The width of the generated gaussian
     * @param amp     The amplitude of the generated gaussian (in %)
     * @param abs     true if you want the absolute values and false otherwise
     */
    private void setShapeTableValuesFromGaussGen(Table table, int nbpoint, double width, double amp, Boolean abs) throws Exception {
        String name = "Gaussian";
        TableGeneratorInterface gen;
        gen = LoaderFactory.getTableGeneratorPluginLoader().getPluginByName(name);
        if (gen == null) {
            throw new IllegalStateException("Table generator not found: " + name);
        }
        table.setGenerator(gen);
        gen.getParams().get(0).setValue(nbpoint);
        gen.getParams().get(1).setValue(width);
        gen.getParams().get(2).setValue(amp);
        gen.getParams().get(3).setValue(abs);//abs

        gen.generate();

    }


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  Offset Fequency
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     * prepare txFrequencyOffsetTable for slice selection
     *
     * @param gradSlice             : Slice selection Gradient
     * @param nbSlice               : number of slice
     * @param spacing_between_slice : space between Slice
     * @param offCenterDistance3D   : offcenter FOV
     */
    public void prepareOffsetFreqMultiSlice(Gradient gradSlice, int nbSlice, double spacing_between_slice, double offCenterDistance3D) {
        numberOfFreqOffset = nbSlice;
        double grad_amp_slice_mTpm = gradSlice.getAmplitude_mTpm();
        double frequencyCenter3D90 = calculateOffsetFreq(grad_amp_slice_mTpm, offCenterDistance3D);
        double slice_thickness = Double.isNaN(gradSlice.getSliceThickness()) ? 0 : gradSlice.getSliceThickness();
        double multi_planar_fov = (numberOfFreqOffset - 1) * (spacing_between_slice + slice_thickness);
        double multiPlanarFreqOffset = multi_planar_fov * grad_amp_slice_mTpm * (GradientMath.GAMMA);
        txFrequencyOffsetTable = new double[numberOfFreqOffset];
        for (int i = 0; i < numberOfFreqOffset; i++) {
            double tx_frequency_offset = (multiPlanarFreqOffset / 2) - (numberOfFreqOffset == 1 ? 0 : i * multiPlanarFreqOffset / (numberOfFreqOffset - 1)) + frequencyCenter3D90;
            txFrequencyOffsetTable[i] = tx_frequency_offset;
        }
    }

    /**
     * reorder txFrequencyOffsetTable acording to the plugin
     *
     * @param acquisitionPointsPerSlice  : acquiered data point in a single scan
     * @param slicesAcquiredInSingleScan : number of acquiered slice in a single scan
     */
    public void reoderOffsetFreq(TransformPlugin plugin, int acquisitionPointsPerSlice, int slicesAcquiredInSingleScan) {
        double sliceNumber;
        double[] offset_table = new double[numberOfFreqOffset];
        for (int k = 0; k < numberOfFreqOffset; k++) {
            int[] indexScan = plugin.invTransf(0, 0, k, 0);
            sliceNumber = indexScan[0] / (acquisitionPointsPerSlice) + indexScan[2] * slicesAcquiredInSingleScan;
            offset_table[(int) sliceNumber] = txFrequencyOffsetTable[k];
        }
        txFrequencyOffsetTable = offset_table;

//
//
//            public int plugin.getSliceScanOrder(int k){
//                return k;
//            }
//            }
//            @Override
//            public int sliceNumber = plugin.getSlicecanOrder(int k){
//                int[] indexScan = invTransf(0, 0, k, 0);
//                sliceNumber = indexScan[0] / (matrix1D * etl) + indexScan[2] * matrix3D / nbshoot;
//                return k
//            }
//


    }

    public double calculateOffsetFreq(double grad_amp_mTpm, double offCenterDistance) {
        return (-grad_amp_mTpm * offCenterDistance * (GradientMath.GAMMA));
    }

    public void setFrequencyOffset(Order order) {
        FrequencyOffsetOrder = order;
        setFrequencyOffset();
    }


    public void setFrequencyOffset(double value) {
        txFrequencyOffsetTable = new double[1];
        txFrequencyOffsetTable[0] = value;
        numberOfFreqOffset = 1;
        setFrequencyOffset();
    }

    public void addFrequencyOffset(double... values) {
        int tmpnumberOfFreqOffset = numberOfFreqOffset;
        if (numberOfFreqOffset != -1) {
            double[] tmpTable = txFrequencyOffsetTable.clone();
            txFrequencyOffsetTable = new double[numberOfFreqOffset + values.length];
            for (int i = 0; i < tmpTable.length; i++) {
                txFrequencyOffsetTable[i] = tmpTable[i];
            }
            numberOfFreqOffset = numberOfFreqOffset + values.length;
        } else {
            numberOfFreqOffset = values.length;
            txFrequencyOffsetTable = new double[numberOfFreqOffset];
        }

        for (int i = 0; i < values.length; i++) {
            txFrequencyOffsetTable[txFrequencyOffsetTable.length - values.length + i] = values[i];
        }
    }

    public void setFrequencyOffset() {
        setSequenceTableValues(FrequencyOffsetTable, FrequencyOffsetOrder);
        if (numberOfFreqOffset != -1) {
            for (int k = 0; k < numberOfFreqOffset; k++) {
                FrequencyOffsetTable.add(txFrequencyOffsetTable[k]);
            }
        } else {
            FrequencyOffsetTable.add(0);
        }
    }

    public void setFrequencyOffset(double... value) {
        numberOfFreqOffset = value.length;
        txFrequencyOffsetTable = new double[numberOfFreqOffset];
        for (int k = 0; k < numberOfFreqOffset; k++) {
            txFrequencyOffsetTable[k] = value[k];
        }
        setFrequencyOffset();
    }

    public void setFrequencyOffset(Order order, double... value) {
        setFrequencyOffset(value);
        setFrequencyOffset(order);
        // System.out.println(txFrequencyOffsetTable[0]);
    }

    public double getFrequencyOffset(int k) {
        return txFrequencyOffsetTable[k];
    }

    /**
     * calculate and set the FrequencyOffset to compensate another pulse
     *
     * @param pulse :  to compensate
     * @param ratio : ratio of the pule duration to be compensated
     */
    public void setCompensationFrequencyOffset(RFPulse pulse, double ratio) {
        FrequencyOffsetOrder = pulse.getFrequencyOffsetOrder();
        numberOfFreqOffset = pulse.getNumberOfFreqOffset();
        if (numberOfFreqOffset != -1) {
            txFrequencyOffsetTable = new double[numberOfFreqOffset];
            for (int k = 0; k < numberOfFreqOffset; k++) {
                txFrequencyOffsetTable[k] = -((pulse.getFrequencyOffset(k) * pulse.getPulseDuration() * ratio) % 1) / pulseDuration;
            }
        }
        setFrequencyOffset();
    }


    /**
     * calculate and set the FrequencyOffset to compensate another pulse
     *
     * @param pulse :  to compensate
     * @param time  : time duration of the pulse to be compensated
     */
    public void setCompensationFrequencyOffsetWithTime(RFPulse pulse, double time) {
        FrequencyOffsetOrder = pulse.getFrequencyOffsetOrder();
        numberOfFreqOffset = pulse.getNumberOfFreqOffset();
        txFrequencyOffsetTable = new double[numberOfFreqOffset];
        for (int k = 0; k < numberOfFreqOffset; k++) {
            txFrequencyOffsetTable[k] = -((pulse.getFrequencyOffset(k) * time) % 1) / this.pulseDuration;
        }
        setFrequencyOffset();
    }

    /**
     * calculate and set the FrequencyOffset for and off center FOV 1D
     *
     * @param grad                :  readout Gradient
     * @param off_center_distance : off_center_distance to be compensated
     */
    public void setFrequencyOffsetReadout(Gradient grad, double off_center_distance) {
        numberOfFreqOffset = 1;
        double grad_amp_read_read_mTpm = grad.getAmplitude_mTpm();// amplitude in T/m
        txFrequencyOffsetTable = new double[numberOfFreqOffset];
        txFrequencyOffsetTable[0] = -grad_amp_read_read_mTpm * off_center_distance * (GradientMath.GAMMA);
        setSequenceTableValues(FrequencyOffsetTable, FrequencyOffsetOrder, txFrequencyOffsetTable[0]);
    }

    public void setFrequencyOffsetReadoutEchoPlanar(Gradient grad, double off_center_distance, int ETL, Order tableorder) {
        numberOfFreqOffset = ETL;
        FrequencyOffsetOrder = tableorder;
        double grad_amp_read_read_mTpm = grad.getAmplitude_mTpm();// amplitude in T/m
        System.out.println("numberOfFreqOffset " + numberOfFreqOffset);
        txFrequencyOffsetTable = new double[numberOfFreqOffset];
        for (int i = 0; i < numberOfFreqOffset; i++) {
            if (i % 2 == 0) {
                txFrequencyOffsetTable[i] = -grad_amp_read_read_mTpm * off_center_distance * (GradientMath.GAMMA);
            } else {
                txFrequencyOffsetTable[i] = grad_amp_read_read_mTpm * off_center_distance * (GradientMath.GAMMA);
            }
        }
        setSequenceTableValues(FrequencyOffsetTable, FrequencyOffsetOrder, txFrequencyOffsetTable);

    }

    /**
     * calculate and set the FrequencyOffset for to induce a phase offset
     *
     * @param angleDegree :  readout Gradient
     */
    public void setFrequencyOffsetForPhaseShift(double angleDegree) {
        numberOfFreqOffset = 1;
        txFrequencyOffsetTable = new double[numberOfFreqOffset];
        txFrequencyOffsetTable[0] = -Math.round((((angleDegree / 360.0) % 1.0) / this.pulseDuration));
        setSequenceTableValues(FrequencyOffsetTable, FrequencyOffsetOrder, txFrequencyOffsetTable[0]);
        return;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  Phase
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public void setPhase(Order order, double... phaseValues) {
        setSequenceTableValues(phase, order, phaseValues);
    }

    public void setPhase(double phaseValues) {
        setSequenceTableSingleValue(phase, phaseValues);
    }

    // ---------------------------------------------------------------
    // ----------------- General Methode----------------------------------------------
    private void setSequenceTableSingleValue(Table table, double... values) {
        // uses Order.One because there are no tables in this dimension: compilation issue
        setSequenceTableValues(table, Order.FourLoop, values);
    }

    private void setSequenceTableValues(Table table, Order order, double... values) {
        table.clear();
        table.setOrder(order);
        table.setLocked(true);
        for (double value : values) {
            table.add(value);
        }
    }

    private double ceilToSubDecimal(double numberToBeRounded, double Order) {
        return Math.ceil(numberToBeRounded * Math.pow(10, Order)) / Math.pow(10, Order);
    }
}
