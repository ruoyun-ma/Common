package common;

import rs2d.commons.log.Log;
import rs2d.spinlab.data.transformPlugin.TransformPlugin;
import rs2d.spinlab.instrument.util.GradientMath;
import rs2d.spinlab.plugins.loaders.LoaderFactory;
import rs2d.spinlab.plugins.loaders.PluginLoaderInterface;
import rs2d.spinlab.sequence.Sequence;
import rs2d.spinlab.sequence.table.Shape;
import rs2d.spinlab.sequence.table.Table;
import rs2d.spinlab.sequence.table.Utility;
import rs2d.spinlab.sequence.table.generator.TableGeneratorInterface;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.param.NumberEnum;
import rs2d.spinlab.tools.param.NumberParam;
import rs2d.spinlab.tools.param.Param;
import rs2d.spinlab.tools.table.Order;
import rs2d.spinlab.tools.utility.Nucleus;
import rs2d.spinlab.api.PowerComputation;
import rs2d.spinlab.api.Hardware;

import java.util.Arrays;
import java.util.List;

import static java.util.stream.IntStream.*;

/**
 * Commons V2.8: add instrument interface
 * Class RFPulse
 * V2.10 correct RF calib
 * V2.9 SLR pulse & constructor V2019.06
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

    private double power_90 = Double.NaN;
    private double power_180 = Double.NaN;
    private double tx_amp = Double.NaN;
    private double flipAngle = Double.NaN;

    private double pulseDuration;
    private double observeFrequency = Double.NaN;
    private Nucleus nucleus = Nucleus.H1;
    private double powerPulse = Double.NaN;

    //XG
    private double sincGenRampSlope = 0.2;
    private Table attOffsetTable = null;//new Table(0, NumberEnum.TxAtt);

    // added parameters for HS pulse
    private double HSbandwidth;
    private double HStruncation;
    private double HSorder;
    private double HSgammaB1;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  Constructor
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public RFPulse() {

    }

    public RFPulse(Param attPar, Table amplitudeTab, Table phaseTab,
                   Table timeTab, Shape shape, Shape shapePhase, Table offsetFreqTab, Nucleus nucleus) {
        amplitudeTable = amplitudeTab;
        attParam = attPar;
        phase = phaseTab;

        timeTable = timeTab;
        this.shape = shape;
        this.shapePhase = shapePhase;

        pulseDuration = timeTable.get(0).doubleValue();
        FrequencyOffsetTable = offsetFreqTab;

        this.nucleus = nucleus;

    }


    public RFPulse(Table timeTab, Table offsetFreqTab, Nucleus nucleus) {
        timeTable = timeTab;

        pulseDuration = timeTab.get(0).doubleValue();
        FrequencyOffsetTable = offsetFreqTab;

        this.nucleus = nucleus;

    }

    public RFPulse(Table timeTab, Table offsetFreqTab, Table phaseTab, Nucleus nucleus) {
        timeTable = timeTab;
        phase = phaseTab;
        pulseDuration = timeTab.get(0).doubleValue();
        FrequencyOffsetTable = offsetFreqTab;

        this.nucleus = nucleus;
    }

    //XG
    public RFPulse(Param attPar, Table attOffsetPar, Table amplitudeTab, Table phaseTab,
                   Table timeTab, Shape shape, Shape shapePhase, Table offsetFreqTab, Nucleus nucleus) {
        //new RFPulse(attPar, amplitudeTab, phaseTab, timeTab, shape, shapePhase, offsetFreqTab);
        amplitudeTable = amplitudeTab;
        attParam = attPar;
        phase = phaseTab;

        timeTable = timeTab;
        this.shape = shape;
        this.shapePhase = shapePhase;

        pulseDuration = timeTable.get(0).doubleValue();
        FrequencyOffsetTable = offsetFreqTab;
        if (attOffsetPar != null) {
            attOffsetTable = attOffsetPar;
        }

        this.nucleus = nucleus;

    }

    public static RFPulse createRFPulse(Sequence sequence, GeneratorSequenceParamEnum txAttParam, GeneratorSequenceParamEnum amplitudeTab, GeneratorSequenceParamEnum txPhaseTab,
                                        GeneratorSequenceParamEnum timeTab, GeneratorSequenceParamEnum shape, GeneratorSequenceParamEnum shapePhaseShape, GeneratorSequenceParamEnum offsetTab, Nucleus nucleus) {
        return new RFPulse(sequence.getPublicParam(txAttParam.name()), sequence.getTable(amplitudeTab.name()), sequence.getTable(txPhaseTab.name()),
                sequence.getPublicTable(timeTab.name()), (Shape) sequence.getTable(shape.name()), (Shape) sequence.getTable(shapePhaseShape.name()), sequence.getTable(offsetTab.name()), nucleus);
    }

    //XG
    public static RFPulse createRFPulse(Sequence sequence, GeneratorSequenceParamEnum txAttParam, GeneratorSequenceParamEnum txAttOffsetParam, GeneratorSequenceParamEnum amplitudeTab, GeneratorSequenceParamEnum txPhaseTab,
                                        GeneratorSequenceParamEnum timeTab, GeneratorSequenceParamEnum shape, GeneratorSequenceParamEnum shapePhaseShape, GeneratorSequenceParamEnum offsetTab, Nucleus nucleus) {
        return new RFPulse(sequence.getPublicParam(txAttParam.name()), sequence.getPublicTable(txAttOffsetParam.name()), sequence.getTable(amplitudeTab.name()), sequence.getTable(txPhaseTab.name()),
                sequence.getPublicTable(timeTab.name()), (Shape) sequence.getTable(shape.name()), (Shape) sequence.getTable(shapePhaseShape.name()), sequence.getTable(offsetTab.name()), nucleus);
    }

    public static RFPulse createRFPulse(Sequence sequence, GeneratorSequenceParamEnum timeTab, GeneratorSequenceParamEnum offsetTab, GeneratorSequenceParamEnum txPhaseTab, Nucleus nucleus) {
        return new RFPulse(sequence.getPublicTable(timeTab.name()), sequence.getTable(offsetTab.name()), sequence.getTable(txPhaseTab.name()), nucleus);
    }

    public static RFPulse createRFPulse(Sequence sequence, GeneratorSequenceParamEnum timeTab, GeneratorSequenceParamEnum offsetTab, Nucleus nucleus) {
        return new RFPulse(sequence.getPublicTable(timeTab.name()), sequence.getTable(offsetTab.name()), nucleus);
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

    public double getFlipAngle() {
        if (Double.isNaN(flipAngle))
            calculateFlipAngle();
        return flipAngle;
    }

    public double getPower() {
        if (Double.isNaN(powerPulse))
            calculatePower();
        return powerPulse;
    }

    public Table getAmplitudeTable() {
        return amplitudeTable;
    }

    public Order getFrequencyOffsetOrder() {
        return FrequencyOffsetOrder;
    }

    public Table getPhaseTable(){ return phase;}

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

    public NumberParam getAttParam() {
        return (NumberParam) attParam;
    }

    public void setAttOffset(int att) {
        setSequenceTableSingleValue(attOffsetTable, att);
    }

    public void setAttOffset(Table attTable) {
        setSequenceTableSingleValue(attOffsetTable, attTable.get(0).intValue());
    }

    public boolean createAttOffset(Sequence sequence, GeneratorSequenceParamEnum txAttOffsetParam) {
        if (sequence.getPublicTable(txAttOffsetParam.name()) != null) {
            attOffsetTable = sequence.getPublicTable(txAttOffsetParam.name());
            setAttOffset(0);
            return true;
        } else
            return false;
    }

    public int getAttOffset() {
        if (attOffsetTable != null)
            return attOffsetTable.get(0).intValue();
        else
            return -1;
    }

    public void setAmp(double amp) {
        tx_amp = amp;
        setSequenceTableSingleValue(amplitudeTable, tx_amp);
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
    }

    public void setHSbandwidth (double bw) {HSbandwidth = bw; }

    public void setHStruncation (double tc) {HStruncation = tc; }

    public void setHSorder (double od) {HSorder = od;}

    public void setHSgammaB1 (double gB1) {HSgammaB1 = gB1;}

    public double getHSbandwidth() {return HSbandwidth;}

    public double getHStruncation() {return HStruncation;}

    public double getHSorder() {return HSorder;}

    public double getHSgammaB1() {return HSgammaB1;}
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  Calculation Amp Att...
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     * Automatic Calibration for a given Flip angle with ATT set to get 180° -> 80%amp.
     * set amplitudeTable and attParam
     *
     * @param flipAngle         : flip angle of the pulse
     * @param observe_frequency :set pulse property
     * @param txRoute           : Tx channel
     * @param nucleus           :set pulse propert
     * @return test_change_time : false if the pulse duration has changed
     */
    public boolean setAutoCalibFor180(double flipAngle, double observe_frequency, List<Integer> txRoute, Nucleus nucleus) {
        this.flipAngle = flipAngle;
        observeFrequency = observe_frequency;
        this.nucleus = nucleus;
        boolean test_change_time = true;
        if (txAtt == -1) {
            test_change_time = prepTxAttFor180(80, txRoute);
        }
        // Calculate amp with the new and real Att
        double tx_amp90 = calculateTxAmp90(txRoute);
        tx_amp = tx_amp90 * flipAngle / 90;
        setSequenceTableSingleValue(amplitudeTable, tx_amp);
        attParam.setValue(txAtt);
        // set calculated parameters to display values & sequence
        return test_change_time;
    }

    /**
     * Automatic Calibration: ATT set to get 180°-> tx_amp%.
     *
     * @param tx_amp  : of the 180° pulse
     * @param txRoute : Tx channel
     * @return test_change_time : false if the pulse duration has changed
     */
    private boolean prepTxAttFor180(double tx_amp, List<Integer> txRoute) {
        boolean test_change_time = true;
        double instrument_length_180 = PowerComputation.getHardPulse180Width(nucleus.name());
        double power_factor = Utility.complexPowerFillingFactor(shape, shapePhase) / getSlrPowerFactors180();       // get RF pulse power factor from instrument to calculate RF pulse amplitude
        double instrument_power_180 = PowerComputation.getHardPulse180Power(nucleus.name()) / power_factor;
        power_180 = instrument_power_180 * Math.pow(instrument_length_180 / pulseDuration, 2);

        if (power_180 > Hardware.getMaxRfPowerPulsed(nucleus.name())) {  // TX LENGTH 90 MIN
            pulseDuration = ceilToSubDecimal(instrument_length_180 / Math.sqrt(Hardware.getMaxRfPowerPulsed(nucleus.name()) / instrument_power_180), 6);
            setSequenceTableSingleValue(timeTable, pulseDuration);
            power_180 = instrument_power_180 * Math.pow(instrument_length_180 / pulseDuration, 2);
            test_change_time = false;
        }
        // Calculate Att to get a 180° RF pulse around 80% amp
        if (txAtt == -1)
            txAtt = PowerComputation.getTxAttenuation(txRoute.get(0), power_180, observeFrequency, tx_amp);
        return test_change_time;
    }

    /**
     * calculate 90° Pulse Amplitude for txAtt
     *
     * @param txRoute : Tx channel
     * @return tx_amp
     */
    private double calculateTxAmp90(List<Integer> txRoute) {
        txAtt = ((NumberParam) attParam).getValue().intValue();
        //XG: adding txAttOffset, 2021/5/4
        if (attOffsetTable != null)
            txAtt += attOffsetTable.get(0).intValue();

        double tx_amp;
        double power_factor = Utility.complexPowerFillingFactor(shape, shapePhase) / getSlrPowerFactors90();       // get RF pulse power factor from instrument to calculate RF pulse amplitude
        double instrument_power_90 = PowerComputation.getHardPulse90Power(nucleus.name()) / power_factor;
        double instrument_length = PowerComputation.getHardPulse90Width(nucleus.name());
        power_90 = instrument_power_90 * Math.pow(instrument_length / pulseDuration, 2);
        tx_amp = PowerComputation.getTxAmplitude(txRoute.get(0), power_90, observeFrequency, txAtt);
        return tx_amp;
    }

    /**
     * calculate 180° Pulse Amplitude for txAtt
     *
     * @param txRoute : Tx channel
     * @return tx_amp
     */
    private double calculateTxAmp180(List<Integer> txRoute) {
        txAtt = ((NumberParam) attParam).getValue().intValue();
        //XG: adding txAttOffset, 2021/5/4
        if (attOffsetTable != null)
            txAtt += attOffsetTable.get(0).intValue();

        double tx_amp;
        double power_factor = Utility.complexPowerFillingFactor(shape, shapePhase) / getSlrPowerFactors180();       // get RF pulse power factor from instrument to calculate RF pulse amplitude
        double instrument_power_180 = PowerComputation.getHardPulse180Power(nucleus.name()) / power_factor;
        double instrument_length = PowerComputation.getHardPulse180Width(nucleus.name());
        power_180 = instrument_power_180 * Math.pow(instrument_length / pulseDuration, 2);
        tx_amp = PowerComputation.getTxAmplitude(txRoute.get(0), power_180, observeFrequency, txAtt);
        return tx_amp;
    }

    /**
     * Calculate the power needed to achieve the flipAngle and check if it exceeds the instrument power limit
     *
     * @param flipAngle        set pulse property
     * @param observeFrequency set pulse property
     * @param nucleus          set pulse property
     * @return test_change_time = false if the pulseDuration was increase because of exceeding power max
     */
    public boolean prepPower(double flipAngle, double observeFrequency, Nucleus nucleus) {
        this.flipAngle = flipAngle;
        this.observeFrequency = observeFrequency; // used to calculate the attenuation and amplitude that match with the power
        this.nucleus = nucleus;
        return calculatePower();
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
        //if (txAtt == -1) { //XG: we not allow to skip powercheck anymore, 2021/3/24
        test_change_time = calculatePower();
        //}
        return test_change_time;
    }

    /**
     * calculate powerPulse from flipAngle and pulseDuration, and check if the power exceed the instrument power limit
     *
     * @return test_change_time = false if the pulseDuration was increase because of exceeding power max
     */
    private boolean calculatePower() {
        boolean b_time_unchanged = true;
        if (Double.isNaN(flipAngle) || Double.isNaN(pulseDuration)) {
            throw new IllegalArgumentException("Power Calculation cannot be done: Flip Angle and/or Pulse Duration are not set!"
                    + "Flip angle is " + flipAngle + ", pulse duration is " + pulseDuration + ".");
        }

        // flip angle < 135° : the power is checked using the 90° hard pulse
        // flip angle > 135° : the power is checked using the 180° hard pulse
        // The RF power is compute using the hard pulse calibration with the closest angle (known relation between power, duration and flip angle) and the shape power factor
        double power_factor = Utility.complexPowerFillingFactor(shape, shapePhase) / (flipAngle < 135 ? getSlrPowerFactors90() : getSlrPowerFactors180());       // get RF pulse power factor from instrument to calculate RF pulse amplitude
        double instrument_length = flipAngle < 135 ? PowerComputation.getHardPulse90Width(nucleus.name()) : PowerComputation.getHardPulse180Width(nucleus.name());
        double instrument_power = (flipAngle < 135 ? PowerComputation.getHardPulse90Power(nucleus.name()) : PowerComputation.getHardPulse180Power(nucleus.name())) / power_factor;
        power_90 = PowerComputation.getHardPulse90Power(nucleus.name()) / power_factor;
        power_180 = PowerComputation.getHardPulse180Power(nucleus.name()) / power_factor;
        powerPulse = instrument_power * Math.pow(instrument_length / pulseDuration, 2) * Math.pow(flipAngle / (flipAngle < 135 ? 90 : 180), 2);

        // If the power exceed the instrument limit increase the pulse duration
        if (powerPulse > Hardware.getMaxRfPowerPulsed(nucleus.name())) {  // TX LENGTH 90 MIN
            pulseDuration = ceilToSubDecimal(instrument_length / Math.sqrt(Hardware.getMaxRfPowerPulsed(nucleus.name()) / (instrument_power * Math.pow(flipAngle / (flipAngle < 135 ? 90 : 180), 2))), 6);
            setSequenceTableSingleValue(timeTable, pulseDuration);
            powerPulse = instrument_power * Math.pow(instrument_length / pulseDuration, 2) * Math.pow(flipAngle / (flipAngle < 135 ? 90 : 180), 2);
            b_time_unchanged = false;
        }

        return b_time_unchanged;
    }

    /**
     *  Calculate the power needed to achieve HSgammaB& and check if it exceeds the instrument power limit
     *
     * @param observe_frequency set pulse property
     * @param nucleus           set pulse property
     * @return test_change_time = false if the pulseDuration was increase because of exceeding power max
     */

    public boolean checkHSPower( double observe_frequency, Nucleus nucleus) {
        observeFrequency = observe_frequency;
        this.nucleus = nucleus;
        boolean test_change_gammaB1 = true;
        test_change_gammaB1 = calculateHSPower();
        return test_change_gammaB1;
    }

    /**
     * calculate HS powerPulse from HSgammaB1, and check if the power exceed the instrument power limit
     *
     * @return test_change_time = false if exceeding power max, HSgammaB1 will be changed
     */
    private boolean calculateHSPower() {
        boolean test_change_gammaB1 = true;

        double TXLength_eff = 0.25/HSgammaB1;
        double instrument_length = PowerComputation.getHardPulse90Width(nucleus.name());
        double instrument_power = PowerComputation.getHardPulse90Power(nucleus.name());
        powerPulse = instrument_power*(instrument_length/TXLength_eff)*(instrument_length/TXLength_eff);
        System.out.println("instrument_length =" + instrument_length) ;
        System.out.println("TXLength_eff =" + TXLength_eff ) ;
        System.out.println("HSgammaB1 =" + HSgammaB1 ) ;
        if (powerPulse > Hardware.getMaxRfPowerPulsed(nucleus.name())) {
            System.out.println("Out of range powerPulse HS = " + powerPulse);
            powerPulse = Hardware.getMaxRfPowerPulsed(nucleus.name());
            TXLength_eff = instrument_length*Math.sqrt(instrument_power/powerPulse);
            HSgammaB1 = 0.25/TXLength_eff;
            test_change_gammaB1 = false;

        }
        return test_change_gammaB1;

    }

    /**
     * Calculate the flip angle obtained with the pulse power knowing the power needed to flip the magnetisation to 90°
     *  (it supposes that the system has been well calibrated)
     * @return flipAngle
     */
    public double calculateFlipAngle() {
        if (Double.isNaN(powerPulse) && !Double.isNaN(powerPulse) && !Double.isNaN(power_90) && !Double.isNaN(pulseDuration)) {
            double instrument_length = PowerComputation.getHardPulse90Width(nucleus.name());
            flipAngle = 90 * Math.sqrt(powerPulse / power_90) * pulseDuration / instrument_length;
        }
        return flipAngle;
    }

    /**
     * calculate and set txAtt in order to get the powerPulse at pulse amplitude = amp
     *
     * @param amp     amplitude for powerPulse
     * @param txRoute Tx channel
     */
    public void prepAtt(double amp, List<Integer> txRoute) {
        // Todo proper handling of exception cases (unreachable Power, No calibration data)
        try {
            txAtt = PowerComputation.getTxAttenuation(txRoute.get(0), powerPulse, observeFrequency, amp);
            txAtt = Math.min(txAtt, 63);

        } catch (Exception e) {
            System.out.println(" No calibration data ");
            Log.info(getClass(), " No calibration data ");
            txAtt = 63;
        }
        attParam.setValue(txAtt);
    }

    /**
     * calculate and set txAtt in order to get the powerPulse at pulse amplitude = amp for a specific power
     *
     * @param refPower  Power used to choose the attenuation
     * @param targetAmp Amplitude for refPower
     * @param txRoute   : Tx channel
     */
    public void prepAtt(double refPower, double targetAmp, List<Integer> txRoute) {

        try {
            txAtt = PowerComputation.getTxAttenuation(txRoute.get(0), refPower, observeFrequency, targetAmp);
            txAtt = Math.min(txAtt, 63);

        } catch (Exception e) {
            System.out.println(" No calibration data ");
            Log.info(getClass(), " No calibration data ");
            txAtt = 63;
        }
        attParam.setValue(txAtt);

    }

    /**
     * prepare and set txAmp according to txAtt , flipAngle
     *
     * @param txRoute : Tx channel
     */
    public void prepTxAmp(List<Integer> txRoute) {
        if (txAtt == -1) {
            txAtt = ((NumberParam) attParam).getValue().intValue();
        }
        tx_amp = PowerComputation.getTxAmplitude(txRoute.get(0), powerPulse, observeFrequency, txAtt);
        setSequenceTableSingleValue(amplitudeTable, tx_amp);
    }

    /**
     * prepare TX amplitude for HS pulse
     *
     *@param txRoute
     * @return tx_amp
     */
    public double prepTxAmpHS(List<Integer> txRoute) {

        txAtt = ((NumberParam) attParam).getValue().intValue();

        double TXLength_eff = 0.25/HSgammaB1;

        double instrument_length = PowerComputation.getHardPulse90Width(nucleus.name());
        double instrument_power = PowerComputation.getHardPulse90Power(nucleus.name());
        double tx_amp_hardpulse = PowerComputation.getTxAmplitude(txRoute.get(0), instrument_power, observeFrequency, txAtt);
        tx_amp = tx_amp_hardpulse*instrument_length/TXLength_eff;
        setSequenceTableSingleValue(amplitudeTable, tx_amp);
        return tx_amp;
    }

    /**
     * prepare and set txAmp according to txAtt , flipAngle
     *
     * @param txRoute : Tx channel
     * @return tx_amp
     */
    public double prepTxAmpMultiFA(List<Integer> txRoute, double[] FA_list, Order order) {
//        //XG: redundancy removed, 2021/5/4
//        //if (txAtt == -1) { //XG: we invoke attParam every time, 2021/3/24
//        txAtt = ((NumberParam) attParam).getValue().intValue();
//        //}
//        //XG: adding txAttOffset, 2021/5/4
//        txAtt += attOffsetTable.get(0).intValue();

        double tx_amp90 = calculateTxAmp90(txRoute);
        double tx_amp180 = calculateTxAmp180(txRoute);
        setSequenceTableValues(amplitudeTable, order);
        range(0, FA_list.length).forEach(i -> {
            tx_amp = (FA_list[i] < 135 ? tx_amp90 : tx_amp180) * FA_list[i] / (FA_list[i] < 135 ? 90 : 180);
            amplitudeTable.add(tx_amp);
        });
        return tx_amp;
    }

    public double prepTxAmpMultiFA(List<Integer> txRoute, List<Double> FA_list, Order order) {
        return prepTxAmpMultiFA(txRoute, FA_list.stream().mapToDouble(d -> d).toArray(), order);
    }

    /**
     * prepare and set txAmp according to GAMMAB1 for adiabatic pulses
     *
     * @param txRoute  : Tx channel
     * @return tx_amp : amplitude
     */
    public double prepTxAmpMultiHS(List<Integer> txRoute, List<Double> gB1_list, Order order) {
        return prepTxAmpMultiHS(txRoute, gB1_list.stream().mapToDouble(d -> d).toArray(), order);
    }
    public double prepTxAmpMultiHS(List<Integer> txRoute, double[] gB1_list, Order order) {

        txAtt = ((NumberParam) attParam).getValue().intValue();

        double instrument_length = PowerComputation.getHardPulse90Width(nucleus.name());
        double instrument_power = PowerComputation.getHardPulse90Power(nucleus.name());
        System.out.println("prepTxAmpHS txATT = " + txAtt);
        double tx_amp_hardpulse = PowerComputation.getTxAmplitude(txRoute.get(0), instrument_power, observeFrequency, txAtt);
        amplitudeTable.clear();
        double TXLength_eff;
        for (int i = 0; i < gB1_list.length; i++){
            if (gB1_list[i] == 0){
                tx_amp = 0.00000;
            } else {
                TXLength_eff = 0.25 / gB1_list[i];
                tx_amp = tx_amp_hardpulse*instrument_length/TXLength_eff;
            }
            amplitudeTable.add(tx_amp);
        } ;

        return tx_amp;
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  Shape
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     * Prepare the RF pulse shape
     *
     * @param pulseName     :GAUSSIAN , SINC3 , SINC5 , HARD
     * @param numberOfPoint The number of point of the generated shape
     * @param type          or "90 degree" "Refocusing (spin-echo)"
     */
    public void setShape(String pulseName, int numberOfPoint, String type) throws Exception {
        shapePhase.clear();
        if ("GAUSSIAN".equalsIgnoreCase(pulseName)) {
            setShapeTableValuesFromGaussGen(shape, numberOfPoint, 0.250, 100, false);
        } else if ("SINC3".equalsIgnoreCase(pulseName)) {
            setShapeTableValuesFromSincGen(shape, numberOfPoint, 2, 100, true, "Hamming");
            for (int i = 0; i < (numberOfPoint); i++) {
                shapePhase.add((i < ((int) (numberOfPoint / 4.0)))
                        || (i >= ((int) (numberOfPoint * 3.0 / 4.0))) ? 180 : 0);
            }
        } else if ("SINC5".equalsIgnoreCase(pulseName)) {
            setShapeTableValuesFromSincGen(shape, numberOfPoint, 3, 100, true, "Hamming");
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
            setTableValuesFromSincGenRamp(shape, numberOfPoint, 3, 100, true, "Hamming", /*0.2*/sincGenRampSlope);
            setTableValuesFromSincGenRampPhase(shapePhase, numberOfPoint, 3, 100, true, "Hamming", /*0.2*/sincGenRampSlope);
        } else if ("SLR_8_5152".equalsIgnoreCase(pulseName)) {
            shape.clear();
            shapePhase.clear();
            String bw = "8.5152";
            setTableValuesFromSLRGen(shape, numberOfPoint, 0, 0, type, true, bw);
            setTableValuesFromSLRPhaseGen(shapePhase, numberOfPoint, 0, 0, type, false, bw);
            isSlr = true;
            slrIndex = 0;
        } else if ("SLR_4_2576".equalsIgnoreCase(pulseName)) {
            shape.clear();
            shapePhase.clear();
            String bw = "4.2576";
            //type="Refocusing (spin-echo)";
            setTableValuesFromSLRGen(shape, numberOfPoint, 0, 0, type, true, bw);
            setTableValuesFromSLRPhaseGen(shapePhase, numberOfPoint, 0, 0, type, false, bw);
            isSlr = true;
            slrIndex = 1;
        } else if ("TONE".equalsIgnoreCase(pulseName)) {
            shape.clear();
            shapePhase.clear();
            setTableValuesFromToneGen(shape, numberOfPoint, 3, 100, false, "Hamming", 2, 0.1);
            setTableValuesFromTonePhaseGen(shapePhase, numberOfPoint, 3, 100, false, "Hamming", 2, 0.1);
        } else if ("HS".equalsIgnoreCase(pulseName)){
            shape.clear();
            shapePhase.clear();
            double beta = 5.289/pulseDuration*2;
            setShapeTableValuesFromHSGen(shape, shapePhase, numberOfPoint,pulseDuration,HSbandwidth, beta );
        } else if ("HSHighOrder".equalsIgnoreCase(pulseName)) {
            shape.clear();
            shapePhase.clear();
            setShapeTableValuesFromHSHighOrdergGen(shape, shapePhase, numberOfPoint, pulseDuration, HSbandwidth, HStruncation, HSorder);
        }
    }

    /**
     * Generate a table of amplitude elements with an SLR pulse generator
     *
     * @param table    The table to be set
     * @param nbpoint  The number of point of the generated SLR pulse
     * @param amp      The amplitude of the generated SLR pulse (in %)
     * @param abs      True if you want the absolute values and false otherwise
     * @param bwstring Bandwidth
     */
    private void setTableValuesFromSLRGen(Table table, int nbpoint, int bw, double amp, String type, boolean abs, String bwstring) throws Exception {
        TableGeneratorInterface gen = null;
        gen = loadTableGenerator("SLR");
        if (gen == null) {
            System.out.println(" no SLR Generator installed");
            Log.error(getClass(), " no SLR Generator installed");
        }
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

    /**
     * Generate a table of phase elements with an SLR pulse generator
     *
     * @param table    The table to be set
     * @param nbpoint  The number of point of the generated SLRphase pulse
     * @param amp      The amplitude of the generated SLRphase pulse (in %)
     * @param abs      True if you want the absolute values and false otherwise
     * @param bwstring Bandwidth
     */
    private void setTableValuesFromSLRPhaseGen(Table table, int nbpoint, int bw, double amp, String type, boolean abs, String bwstring) throws Exception {
        TableGeneratorInterface gen = null;
        gen = loadTableGenerator("SLRPhase");
        if (gen == null) {
            System.out.println(" no SLR Generator installed");
            Log.error(getClass(), " no SLR Generator installed");
        }
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

    /**
     * Generate a table of amplitude elements with a Sinus Cardinal ramp generator
     *
     * @param table   The table to be set
     * @param nbpoint The number of point of the generated sinus cardinal ramp
     * @param nblobe  The number of lob of the generated sinus cardinal ramp
     * @param amp     The amplitude of the generated sinus cardinal ramp (in %)
     * @param abs     true if you want the absolute values and false otherwise
     * @param slope   Ramp slope
     */
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

    /**
     * Generate a table of phase elements with a Sinus Cardinal ramp generator
     *
     * @param table   The table to be set
     * @param nbpoint The number of point of the generated sinus cardinal ramp
     * @param nblobe  The number of lob of the generated sinus cardinal ramp
     * @param amp     The amplitude of the generated sinus cardinal ramp (in %)
     * @param abs     true if you want the absolute values and false otherwise
     * @param slope   Ramp slope
     */
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

    private void setTableValuesFromToneGen(Table table, int nbpoint, int nblobe, double amp, Boolean abs, String window, double slope, double slcborder) throws Exception {
        TableGeneratorInterface gen = null;
        gen = loadTableGenerator("Tone");
        gen.getParams().get(0).setValue(nbpoint);
        gen.getParams().get(1).setValue(nblobe);
        gen.getParams().get(2).setValue(amp);
        gen.getParams().get(3).setValue(abs);//abs
        gen.getParams().get(4).setValue(window);//abs
        gen.getParams().get(5).setValue(slope);
        gen.getParams().get(6).setValue(slcborder);

        table.setGenerator(gen);
        if (gen == null) {
            table.clear();
            table.setFirst(100);
        } else {
            gen.generate();
        }
    }


    private void setTableValuesFromTonePhaseGen(Table table, int nbpoint, int nblobe, double amp, Boolean abs, String window, double slope, double slcborder) throws Exception {
        TableGeneratorInterface gen = null;
        gen = loadTableGenerator("TonePhase");
        gen.getParams().get(0).setValue(nbpoint);
        gen.getParams().get(1).setValue(nblobe);
        gen.getParams().get(2).setValue(amp);
        gen.getParams().get(3).setValue(abs);//abs
        gen.getParams().get(4).setValue(window);//abs
        gen.getParams().get(5).setValue(slope);
        gen.getParams().get(6).setValue(slcborder);

        table.setGenerator(gen);
        if (gen == null) {
            table.clear();
            table.setFirst(100);
        } else {
            gen.generate();
        }
    }

    /**
     * Table generator
     *
     * @param generatorName Name of the shape to generate
     * @return Table generator
     */
    private TableGeneratorInterface loadTableGenerator(String generatorName) throws Exception {
        TableGeneratorInterface gen = null;
        PluginLoaderInterface loader = LoaderFactory.getTableGeneratorPluginLoader();
        if (loader.containsPlugin(generatorName)) {
            gen = (TableGeneratorInterface) loader.getPluginByName(generatorName);
        }
        return gen;
    }

    /**
     * Generate a table of elements with a Sinus Cardinal generator
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

    public void setSincGenRampSlope(double val) {
        sincGenRampSlope = val;
    }

    /**
     * Generate a table of element of HyperbolicSecant pulse
     *
     * @param table_shape   The table to be set
     * @param table_shape_phase   The table to be set
     * @param nbpoint The number of point of the generated Hyperbolic secant signal, should be a multiple of 4
     * @param duration Duration of pulse in seconds.
     * @param bandwidthHS bandwidth of the HS pulse (Hz), currently set to 8400 Hz
     * @param beta AM waveform parameter (rad/s)
     *
     *
     */
    private void setShapeTableValuesFromHSGen(Table table_shape, Table table_shape_phase, int nbpoint, double duration, double bandwidthHS, double beta) {
       /* r"""Design a hyperbolic secant adiabatic pulse.

        mu * beta becomes the amplitude of the frequency sweep

        Args:
        nbpoint (int): number of samples (should be a multiple of 4).
                beta (float): AM waveform parameter.
                mu (float): a constant, determines amplitude of frequency sweep.
        duration(float): pulse time (s).
        beta (float): AM waveform parameter. (rad/s)
        bandwidthHS : bandwidth of HS pulse. (Hz)

        Calculates

        - **table_shape** (*array*): AM waveform.
        - **table_shape_phase** (*array*): phase

        References:
        Bernstein M, King K, Zhou X, Handbook of PulseSequence Chapter 6, 6.2.3, pp194
        """
        */

        double mu = bandwidthHS * Math.PI / beta; // mu determines amplitude of frequency sweep
        double[] t = new double[nbpoint];
        //double [] shape_mag = new double[nbpoint];
        double [] shape_pha = new double[nbpoint];
        double [] shape_pha_wrap = new double[nbpoint];
        for (int i = 0; i < nbpoint; i++) {
            double ii;
            ii = i;
            t[i] = (ii - nbpoint / 2) / nbpoint * duration;

            table_shape.add (100 / Math.cosh(t[i] * beta));
            shape_pha[i] = -mu * Math.log(Math.cosh(beta * t[i]));
            shape_pha_wrap[i] = -Math.atan2(Math.sin(shape_pha[i]),Math.cos(shape_pha[i]));
            if (shape_pha_wrap[i] < 0){
                shape_pha_wrap[i] = shape_pha_wrap[i] / Math.PI * 180.0 + 360.0;
            } else {
                shape_pha_wrap[i] = shape_pha_wrap[i] / Math.PI * 180.0;
            }
            table_shape_phase.add(shape_pha_wrap[i]);

        }

    }



    /**
     * Generate a table of element for High Order HS pulse
     *
     * @param table_shape   The table to be set
     * @param table_shape_phase   The table to be set
     * @param nbpoint The number of point of the generated Hyperbolic secant signal, should be a multiple of 4
     * @param duration Duration of pulse in seconds.
     * @param bandwidthHS bandwidth of the HS pulse (Hz), currently set to 8400 Hz
     * @param powerNumber power of HS pulse
     *
     *
     */
    private void setShapeTableValuesFromHSHighOrdergGen(Table table_shape, Table table_shape_phase, int nbpoint, double duration, double bandwidthHS, double tc, double powerNumber) {
       /* r"""Design a hyperbolic secant adiabatic pulse.
       code implementation based on rf design tool from FID_A
       https://github.com/CIC-methods/FID-A/blob/master/rfPulseTools/rf_hs.m

        mu * beta becomes the amplitude of the frequency sweep

        Args:
        nbpoint (int): number of samples (should be a multiple of 4).
                beta (float): AM waveform parameter.
                mu (float): a constant, determines amplitude of frequency sweep.
        duration(float): pulse time (s).

        bandwidthHS : bandwidth of HS pulse. (Hz)
        truncation factor: truncation factor for HS pulses, should be smaller than1

        Calculates

        - **table_shape** (*array*): AM waveform.
        - **table_shape_phase** (*array*): phase

        References:
        Bernstein M, King K, Zhou X, Handbook of PulseSequence Chapter 6, 6.2.3, pp194
        """
        */

        double B = Math.log(1/tc + Math.sqrt(1/tc/tc -1)); // arcsech (tc)
        double A = bandwidthHS / 2;

        double [] tau = new double[nbpoint];
        double [] F1 = new double[nbpoint];
        double [] F1cumsum = new double [nbpoint/2];
        double [] F1cumsumFlip = new double [nbpoint/2];
        double [] F2 = new double[nbpoint];
        double nbpoint_db = nbpoint;

        System.out.println(" -------------HS print start --------------");
        System.out.println("order = "+powerNumber);
        System.out.println("truncation = "+tc);
        System.out.println("bandwidth = "+bandwidthHS);
        System.out.println("B = "+B);
        for (int i = 0; i < nbpoint; i++) {

            double ii;
            ii = i;
            tau[i] = (ii - nbpoint_db / 2) / nbpoint_db * 2;
            F1[i] = 1/Math.cosh(B * Math.pow(tau[i],powerNumber));
            table_shape.add(F1[i]*100.0);
            //System.out.println(Math.pow(tau[i],powerNumber));

            //table_shape.add(100 / Math.cosh(t[i] * beta));
            //table_shape_phase.add( 180 / Math.PI * mu * Math.cosh(beta * t[i]));
        }
        F1cumsum[0] = F1[nbpoint/2];
        for (int i = nbpoint/2+1; i < nbpoint; i++){
            F1cumsum[i-nbpoint/2] = F1cumsum[i-nbpoint/2-1] + F1[i];
            F1cumsumFlip[nbpoint-i-1] = -F1cumsum[i-nbpoint/2];
        }
        for (int i = 0; i < nbpoint/2; i++){
            F2[i+nbpoint/2] = F1cumsum[i];
            F2[i] = F1cumsumFlip[i];
        }
        double F2max = 0;
        for (int i = 0; i < nbpoint; i++){
            if (F2[i] > F2max) {
                F2max = F2[i];
            }
        }
        double dt = duration/nbpoint_db;
        double [] FMcumsum = new double [nbpoint];
        FMcumsum[0] = A * F2[0]/F2max;
        for (int i = 1; i < nbpoint; i++){
            FMcumsum[i] = A * F2[i]/F2max + FMcumsum[i-1];
        }
        double [] ph = new double [nbpoint];
        for (int i = 0; i < nbpoint; i++){
            table_shape_phase.add( FMcumsum[i] *dt *360.0);
        }
    }

    public void fillShape(List<Number> rfShape, List<Number> rfShapePhase){
        this.shape.clear();
        this.shapePhase.clear();
        for (int i = 0; i<rfShape.size(); i++) {
            this.shape.add(rfShape.get(i));
            this.shapePhase.add(rfShapePhase.get(i));
        }
    }

    public void fillShape(List<Number> rfShape){
        this.shape.clear();
        for (int i = 0; i<rfShape.size(); i++) {
            this.shape.add(rfShape.get(i));

        }
        this.shapePhase.clear();
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  Offset Frequency
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
        double multiPlanarFreqOffset = multi_planar_fov * grad_amp_slice_mTpm * (GradientMath.GAMMA / nucleus.getRatio());
        txFrequencyOffsetTable = new double[numberOfFreqOffset];
        for (int i = 0; i < numberOfFreqOffset; i++) {
            double tx_frequency_offset = (multiPlanarFreqOffset / 2) - (numberOfFreqOffset == 1 ? 0 : i * multiPlanarFreqOffset / (numberOfFreqOffset - 1)) + frequencyCenter3D90;
            txFrequencyOffsetTable[i] = tx_frequency_offset;
        }
    }

    /**
     * reorder txFrequencyOffsetTable according to the plugin for interleaved multi-slices acquisition
     *
     * @param acquisitionPointsPerSlice  : acquired data point in a single scan
     * @param slicesAcquiredInSingleScan : number of acquired slice in a single scan
     */
    public void reoderOffsetFreq(TransformPlugin plugin, int acquisitionPointsPerSlice, int slicesAcquiredInSingleScan) {
        double sliceNumber;
        double[] offset_table = new double[numberOfFreqOffset];
        for (int k = 0; k < numberOfFreqOffset; k++) {
            int[] indexScan = plugin.invTransf(0, 0, k, 0);
            // slicesAcquiredInSingleScan are acquired in one TR so length index[2] = number of TR need to acquired all the slices
            // For the nth TR : index[0] = nTR * slicesAcquiredInSingleScan
            sliceNumber = (int) (indexScan[0] / (double) acquisitionPointsPerSlice) + indexScan[2] * slicesAcquiredInSingleScan;
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

    public void reoderOffsetFreq(int slicesAcquiredInSingleScan) {
        if (Math.floor(numberOfFreqOffset / slicesAcquiredInSingleScan) * slicesAcquiredInSingleScan != numberOfFreqOffset)
            Log.error(getClass(), "floor(numberOfFreqOffset / slicesAcquiredInSingleScan) * slicesAcquiredInSingleScan != numberOfFreqOffset");

        double[] offset_table = new double[numberOfFreqOffset];
        int[] sliceNumber_table = new int[numberOfFreqOffset];
        for (int i = 0; i < Math.floor(numberOfFreqOffset / slicesAcquiredInSingleScan); i++) {
            for (int j = 0; j < slicesAcquiredInSingleScan; j++) {
                sliceNumber_table[slicesAcquiredInSingleScan * i + j] = i + (int) Math.floor(numberOfFreqOffset / slicesAcquiredInSingleScan) * j;
            }
        }

        for (int k = 0; k < numberOfFreqOffset; k++) {
            offset_table[k] = txFrequencyOffsetTable[sliceNumber_table[k]];
        }
        txFrequencyOffsetTable = offset_table;
    }

    /**
     * Calculate the frequency offset that corresponds to the off-center distance wanted given the gradient amplitude
     *
     * @param grad_amp_mTpm     : Amplitude of the gradient [mT/m]
     * @param offCenterDistance : Distance to the iso-center of the gradient
     * @return frequency offset
     */

    public double calculateOffsetFreq(double grad_amp_mTpm, double offCenterDistance) {
        return (-grad_amp_mTpm * offCenterDistance * (GradientMath.GAMMA / nucleus.getRatio()));
    }

    /**
     * Choose the loop that will increment the frequency offset and prepare FrequencyOffsetTable with the existing txFrequencyOffsetTable values
     *
     * @param order : loop order
     */

    public void setFrequencyOffset(Order order) {
        FrequencyOffsetOrder = order;
        setFrequencyOffset();
    }

    /**
     * Set txFrequencyOffsetTable to a specific value, and then prepare the FrequencyOffsetTable
     *
     * @param value : value of the frequency offset
     */
    public void setFrequencyOffset(double value) {
        txFrequencyOffsetTable = new double[1];
        txFrequencyOffsetTable[0] = value;
        numberOfFreqOffset = 1;
        setFrequencyOffset();
    }

    /**
     * Add frequency offsets at the end of tx offset array
     *
     * @param values : array of additional frequency offsets
     */
    public void addFrequencyOffset(double... values) {
        if (numberOfFreqOffset != -1) {
            double[] tmpTable = txFrequencyOffsetTable.clone();
            txFrequencyOffsetTable = new double[numberOfFreqOffset + values.length];
            System.arraycopy(tmpTable, 0, txFrequencyOffsetTable, 0, tmpTable.length);
            numberOfFreqOffset = numberOfFreqOffset + values.length;
        } else {
            numberOfFreqOffset = values.length;
            txFrequencyOffsetTable = new double[numberOfFreqOffset];
        }
        System.arraycopy(values, 0, txFrequencyOffsetTable, txFrequencyOffsetTable.length - values.length, values.length);
    }

    /**
     * Copy the values of txFrequencyOffsetTable into FrequencyOffsetTable and set the incremental loop order to FrequencyOffsetOrder
     */
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

    /**
     * Set txFrequencyOffsetTable, and then prepare the FrequencyOffsetTable values
     *
     * @param value : array of frequency offsets
     */
    public void setFrequencyOffset(double... value) {
        numberOfFreqOffset = value.length;
        txFrequencyOffsetTable = new double[numberOfFreqOffset];
        System.arraycopy(value, 0, txFrequencyOffsetTable, 0, numberOfFreqOffset);
        setFrequencyOffset();
    }

    /**
     * Prepare the FrequencyOffsetTable values and loop order
     *
     * @param order : loop order
     * @param value : array of frequency offsets
     */
    public void setFrequencyOffset(Order order, double... value) {
        setFrequencyOffset(value);
        setFrequencyOffset(order);
    }

    /**
     * Return the element at a given index from the txFrequencyOffsetTable Array
     *
     * @param k : array index
     * @return : k-th value of the txFrequencyOffsetTable
     */
    public double getFrequencyOffset(int k) {
        return txFrequencyOffsetTable[k];
    }

    /**
     * Calculate and set the FrequencyOffset to compensate another pulse
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
        txFrequencyOffsetTable[0] = -grad_amp_read_read_mTpm * off_center_distance * (GradientMath.GAMMA / nucleus.getRatio());
        setSequenceTableValues(FrequencyOffsetTable, FrequencyOffsetOrder, txFrequencyOffsetTable[0]);
    }

    /**
     * Prepare the FrequencyOffsetTable for EPI readout (gradient inverted one echo over 2)
     *
     * @param grad                : gradient amplitude
     * @param off_center_distance : Distance to the iso-center of the gradient to be compensated
     * @param ETL                 : Echo train length (number of echoes that compose the echo train)
     * @param tableorder          : Loop order
     */
    public void setFrequencyOffsetReadoutEchoPlanar(Gradient grad, double off_center_distance, int ETL, Order tableorder) {
        numberOfFreqOffset = ETL;
        FrequencyOffsetOrder = tableorder;
        double grad_amp_read_read_mTpm = grad.getAmplitude_mTpm();// amplitude in T/m
        txFrequencyOffsetTable = new double[numberOfFreqOffset];
        for (int i = 0; i < numberOfFreqOffset; i++) {
            if (i % 2 == 0) {
                txFrequencyOffsetTable[i] = -grad_amp_read_read_mTpm * off_center_distance * (GradientMath.GAMMA / nucleus.getRatio());
            } else {
                txFrequencyOffsetTable[i] = grad_amp_read_read_mTpm * off_center_distance * (GradientMath.GAMMA / nucleus.getRatio());
            }
        }
        setSequenceTableValues(FrequencyOffsetTable, FrequencyOffsetOrder, txFrequencyOffsetTable);

    }

    /**
     * calculate and set the FrequencyOffset to induce a phase offset
     *
     * @param angleDegree : readout Gradient
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

    /**
     * Prepare the phase values and loop order
     *
     * @param order       : loop order
     * @param phaseValues : array of phase values
     */
    public void setPhase(Order order, double... phaseValues) {
        setSequenceTableValues(phase, order, phaseValues);
    }

    /**
     * Prepare the phase values, the fourth loop will be used
     *
     * @param phaseValues : array of phase values
     */
    public void setPhase(double phaseValues) {
        setSequenceTableSingleValue(phase, phaseValues);
    }

    // -------------------------------------------------------------------------------
    // ----------------- General Methode----------------------------------------------

    /**
     * Configure the sequence parameter table with specified values, the loop order is set to FourLoop
     *
     * @param table  : Table name
     * @param values : Table values
     */
    private void setSequenceTableSingleValue(Table table, double... values) {
        // uses Order.One because there are no tables in this dimension: compilation issue
        setSequenceTableValues(table, Order.FourLoop, values);
    }

    /**
     * Configure the sequence parameter table with specified values and loop order
     *
     * @param table  : Table name
     * @param order  : Loop order
     * @param values : Table values
     */
    private void setSequenceTableValues(Table table, Order order, double... values) {
        table.clear();
        table.setOrder(order);
        table.setLocked(true);
        for (double value : values) {
            table.add(value);
        }
    }

    /**
     * Ceil a number to a given decimal places
     *
     * @param numberToBeRounded : Double number
     * @param Order             : Digits kept after the decimal point
     */
    private double ceilToSubDecimal(double numberToBeRounded, double Order) {
        return Math.ceil(numberToBeRounded * Math.pow(10, Order)) / Math.pow(10, Order);
    }

    public int getAttParamValue() {
        return ((NumberParam) attParam).getValue().intValue();
    }
}
