package rs2d.sequence.common;

import rs2d.spinlab.data.transformPlugin.TransformPlugin;
import rs2d.spinlab.hardware.controller.HardwareHandler;
import rs2d.spinlab.instrument.util.GradientMath;
import rs2d.spinlab.sequence.Sequence;
import rs2d.spinlab.sequence.table.Shape;
import rs2d.spinlab.sequence.table.Table;
import rs2d.spinlab.sequence.table.Utility;
import rs2d.spinlab.tools.table.Order;

/**
 * Class Gradient
 * V2.1- 2017-10-24 JR
 */
public class Gradient {
    private Table amplitudeTable = null;
    private Table flatTimeTable = null;
    private Shape shapeUpTable = null;
    private Shape shapeDownTable = null;
    private Table rampTimeUpTable;
    private Table rampTimeDownTable;

    private double amplitude = Double.NaN;
    private double staticArea = Double.NaN;

    private double[] amplitudeArray = null;
    private double maxAreaPE = Double.NaN;

    private int steps = -1;
    private Order order = Order.FourLoop;

    private double grad_shape_rise_time = Double.NaN;
    private double equivalentTime = Double.NaN;

    // Slc
    private double sliceThicknessExcitation = Double.NaN;
    private double txBandwidth = Double.NaN;

    // RO
    private double spectralWidth = Double.NaN;

    // PE
    private boolean bPhaseEncoding = false;
    private double fovPhase = Double.NaN;
    private boolean isKSCentred = false;
    private double spoilerExcess = Double.NaN;

    private boolean bRefocalizeGradient = false;


    private static double gMax = GradientMath.getMaxGradientStrength();

    public Gradient(Table amplitudeTab, Table flat_TimeTab, Shape shapeUpTab, Shape shapeDownTab, Table rampTimeUpTab, Table rampTimeDownTab) {
        amplitudeTable = amplitudeTab;
        flatTimeTable = flat_TimeTab;
        shapeUpTable = shapeUpTab;
        shapeDownTable = shapeDownTab;
        rampTimeUpTable = rampTimeUpTab;
        rampTimeDownTable = rampTimeDownTab;
        gMax = GradientMath.getMaxGradientStrength();
        computeShapeRiseTime();
        prepareEquivalentTime();
    }

    public static Gradient createGradient(Sequence sequence, String amplitudeTab, String flat_TimeTab, String shapeUpTab, String shapeDownTab, String rampTimeTab) {
        return new Gradient(sequence.getPublicTable(amplitudeTab), sequence.getPublicTable(flat_TimeTab), (Shape) sequence.getPublicTable(shapeUpTab),
                (Shape) sequence.getPublicTable(shapeDownTab), sequence.getPublicTable(rampTimeTab), sequence.getPublicTable(rampTimeTab));
    }

    public static Gradient createGradient(Sequence sequence, String amplitudeTab, String flat_TimeTab, String shapeUpTab, String shapeDownTab, String rampTimeUpTab, String rampTimeDownTab) {
        return new Gradient(sequence.getPublicTable(amplitudeTab), sequence.getPublicTable(flat_TimeTab), (Shape) sequence.getPublicTable(shapeUpTab),
                (Shape) sequence.getPublicTable(shapeDownTab), sequence.getPublicTable(rampTimeUpTab), sequence.getPublicTable(rampTimeDownTab));
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  general  methodes
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    public double getStaticArea() {
        if (Double.isNaN(staticArea)) {
            return 0.0;
        } else {
            return staticArea;
        }
    }

    public double getAmplitude() {
        if (Double.isNaN(amplitude)) {
            return 0.0;
        } else {
            return amplitude;
        }
    }

    public double getAmplitude_mTpm() {
        if (Double.isNaN(amplitude)) {
            return 0.0;
        } else {
            return amplitude * gMax / 100.0;
        }
    }

    public double getAmplitudeArray(int pos) {
        if (amplitudeArray == null || pos >= steps) {
            return Double.NaN;
        } else {
            return amplitudeArray[pos];
        }
    }

    public double getSliceThickness() {
        return sliceThicknessExcitation;
    }

    public double getSpoilerExcess() {
        return spoilerExcess;
    }


    public int getSteps() {
        return steps;
    }

    public Order getOrder() {
        return order;
    }

    public double getEquivalentTime() {
        return equivalentTime;
    }

    public double getGradShapeRiseTime() {
        computeShapeRiseTime();
        return equivalentTime;
    }

    public double getTotalArea() {
        double totalArea = 0.0;
        if (!Double.isNaN(maxAreaPE)) {
            totalArea += maxAreaPE;
        }
        if (!Double.isNaN(staticArea)) {
            totalArea += staticArea;
        }
        return totalArea;
    }

    public double calculateStaticArea() {
        if (Double.isNaN(equivalentTime)) {
            prepareEquivalentTime();
        }
        staticArea = (equivalentTime) * amplitude;
        return staticArea;
    }

    /**
     * calculate Static Amplitude from staticArea
     */
    public double calculateStaticAmplitude() {
        if (Double.isNaN(equivalentTime)) {
            prepareEquivalentTime();
        }
        amplitude = staticArea / (equivalentTime);
        return amplitude;
    }

    /**
     * write the prepared amplitude into the Sequence_parameter
     */
    public void applyAmplitude() {
        double offset = 0.0;
        if (!Double.isNaN(amplitude)) {
            offset = amplitude;
        }
        if (amplitudeArray != null) {
            setSequenceTableValues(amplitudeTable, order);
            for (int i = 0; i < steps; i++) {
                amplitudeTable.add(amplitudeArray[i] + offset);
            }
        } else {
            order = Order.FourLoop;
            setSequenceTableValues(amplitudeTable, order);
            amplitudeTable.add(offset);
        }
    }

    public void applyAmplitude(Order taborder) {
        order = taborder;
        applyAmplitude();
    }

    /**
     * Check if the statid and dynamic gradient do not exceed the gradient max:
     */
    public double[] checkGradientMax() {
        double gradMax = 0;
        double gradMin = 0;
        double offset = 0.0;
        if (!Double.isNaN(amplitude)) {
            offset = amplitude;
        }
        if (amplitudeArray != null) {
            for (int i = 0; i < steps; i++) {
                gradMax = Math.max(gradMax, (amplitudeArray[i] + offset));
                gradMin = Math.min(gradMin, (amplitudeArray[i] + offset));
            }
        } else {
            gradMax = Math.max(gradMax, Math.abs(+offset));
        }
        if ((gradMax > 100.0) || (gradMin < -100.0)) {
            System.out.println("Warning - Gradient Min and  Max: " + gradMax + " " + gradMin);
        }
        double[] gradMaxMin = new double[2];
        gradMaxMin[0] = gradMax;
        gradMaxMin[1] = gradMin;
        return (gradMaxMin);
    }

    /**
     * Calculate equivalent shape rise time from rampTime ans shapes
     * equivalent time for rectangular gradient with same Area and amplitude
     *
     * @return grad_shape_rise_time :
     */
    private double computeShapeRiseTime() {
        double grad_shape_rise_factor_up = Utility.voltageFillingFactor(shapeUpTable);
        double grad_shape_rise_factor_down = Utility.voltageFillingFactor(shapeDownTable);
        double grad_up_rise_time = rampTimeUpTable.get(0).doubleValue();
        double grad_down_rise_time = rampTimeDownTable.get(0).doubleValue();
        grad_shape_rise_time = grad_shape_rise_factor_up * grad_up_rise_time + grad_shape_rise_factor_down * grad_down_rise_time;
        return grad_shape_rise_time;
    }

    /**
     * Calculate equivalentTime of a rectangular gradient with same Area and amplitude
     *
     * @return equivalentTime :
     */
    public double prepareEquivalentTime() {
        if (grad_shape_rise_time == Double.NaN) {
            computeShapeRiseTime();
        }
        equivalentTime = (flatTimeTable.get(0).doubleValue() + grad_shape_rise_time);
        return equivalentTime;
    }

    public void refocalizeGradient() {
        calculateStaticAmplitude();
    }

    public void refocalizeGradient(Gradient grad, double ratio) {
        bRefocalizeGradient = true;
        staticArea = -grad.getStaticArea() * ratio;
        calculateStaticAmplitude();
    }

    public void rePrepare() {
        prepareEquivalentTime();
        if (bPhaseEncoding)
            preparePhaseEncoding();
        if (bRefocalizeGradient)
            refocalizeGradient();

    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //     RO             Readout  methodes             RO
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     * calculate readout Gradient Amplitude
     *
     * @param spectralWidth
     * @param fov
     * @return testSpectralWidth : false if the Spectralwidth need to be nicreased(call getSpectralWidth() )
     */
    public boolean calculateReadoutGradient(double spectralWidth, double fov) throws Exception {
        boolean testSpectralWidth = true;
        this.spectralWidth = spectralWidth;
        amplitude = spectralWidth / ((GradientMath.GAMMA) * fov) * 100.0 / gMax;                 // amplitude in T/m
        if (amplitude > 100.0) {
            this.spectralWidth = solveSpectralWidthMax(fov);
            amplitude = this.spectralWidth / ((GradientMath.GAMMA) * fov) * 100.0 / gMax;                 // amplitude in T/m
            testSpectralWidth = false;
        }
        calculateStaticArea();
        return testSpectralWidth;
    }

    /**
     * calculate spectral width max for a given FOV
     *
     * @param fov
     * @return sw
     */
    public double solveSpectralWidthMax(double fov) throws Exception {
        double spectralWidth = getInferiorSpectralWidth(gMax * GradientMath.GAMMA * fov);
        spectralWidth = 3906250.0 / (Math.ceil(3906250.0 / spectralWidth) + 1); // to get the nearest SW bellow the limit
        spectralWidth = HardwareHandler.getInstance().getSequenceHandler().getCompiler().getNearestSW(spectralWidth);
        return spectralWidth;
    }

    public double getSpectralWidth() {
        return spectralWidth;
    }

    /**
     * set READOUT gradient Amplitude ETL values with back and forth(+/-) sign
     *
     * @param ETL        : number of gradient values
     * @param tableorder
     */
    public void applyReadoutEchoPlanarAmplitude(int ETL, Order tableorder) {
        order = tableorder;
        if (!Double.isNaN(amplitude)) {
            steps = ETL;
            amplitudeArray = new double[steps];
            for (int i = 0; i < steps; i++) {
                if (i % 2 == 0) {
                    amplitudeArray[i] = amplitude;
                } else {
                    amplitudeArray[i] = -amplitude;
                }
            }
            amplitude = 0.0;
        }
        applyAmplitude();
    }

    /**
     * calculate READOUT refocusing gradient Amplitude handeling ETL
     *
     * @param grad : Readout Gradient
     */
    public void refocalizeReadoutGradient(Gradient grad, double ratio) {
        int rOSteps = grad.getSteps();
        if (rOSteps > 0) {
            ratio = (rOSteps % 2) == 1 ? ratio : 1 - ratio;
            staticArea = -grad.getAmplitudeArray(rOSteps - 1) * grad.getEquivalentTime() * ratio;
        } else {
            staticArea = -grad.getAmplitude() * grad.getEquivalentTime() * ratio;
        }
        calculateStaticAmplitude();
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //     Slc             Slice Selection             Slc
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public boolean prepareSliceSelection(double tx_bandwidth, double slice_thickness_excitation) {
        boolean testSliceThickness = true;
        txBandwidth = tx_bandwidth;
        this.sliceThicknessExcitation = slice_thickness_excitation;
        amplitude = (tx_bandwidth / ((GradientMath.GAMMA) * sliceThicknessExcitation)) * 100.0 / gMax;                 // amplitude in T/m
        if (amplitude > 100.0) {
            sliceThicknessExcitation = (tx_bandwidth / ((GradientMath.GAMMA) * gMax));
            amplitude = 100;
            testSliceThickness = false;
        }
        calculateStaticArea();
        return testSliceThickness;
    }

    public void applyAmplitude(double sliceThickness) {
        prepareSliceSelection(txBandwidth, sliceThickness);
        applyAmplitude();
    }


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //     PE             Phase Encoding             PE
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    public void preparePhaseEncoding() {
        preparePhaseEncoding(steps, fovPhase, isKSCentred);
    }

    public void preparePhaseEncoding(int matrixDimension, double fovDimPE, boolean isKSCentred) {
        bPhaseEncoding = true;
        steps = matrixDimension;
        fovPhase = fovDimPE;
        this.isKSCentred = isKSCentred;
        double grad_total_area_phase = prepPhaseGradTotalArea(steps, fovPhase);
        double grad_index_max_phase = prepPhaseGradIndexMax(this.isKSCentred);
        double grad_total_amp_phase = grad_total_area_phase / equivalentTime;
        amplitudeArray = new double[steps];
        for (int i = 0; i < steps; i++) {
            amplitudeArray[i] = -(grad_index_max_phase * grad_total_amp_phase) + i * grad_total_amp_phase / (steps - 1);
        }
    }

    public void preparePhaseEncodingForCheck(int matrixDimensionForCheck, int matrixDimension, double fovDim, boolean isKSCentred) {
        double grad_total_area_phase = prepPhaseGradTotalArea(matrixDimensionForCheck, fovDim);
        double grad_index_max_phase = prepPhaseGradIndexMax(isKSCentred);
        maxAreaPE = grad_index_max_phase * grad_total_area_phase;

        preparePhaseEncoding(matrixDimension, fovDim, isKSCentred);
    }

    public double prepPhaseGradTotalArea(int matrixDimension, double fovPhase) {
        return ((matrixDimension - 1) / ((GradientMath.GAMMA) * fovPhase)) * 100.0 / gMax;
    }

    public double prepPhaseGradIndexMax(boolean isKSCentred) {
        double gradIndexMaxPhase;
        if (isKSCentred) {
            gradIndexMaxPhase = 1 / 2.0;// symetric k'space around zero
        } else {
            gradIndexMaxPhase = 1 / 2.0 + ((steps + 1) % 2) / (2.0 * ((float) steps - 1));// always go trough k0
        }
        return gradIndexMaxPhase;
    }


    public double[] refocalizePhaseEncodingGradient(Gradient grad) {
        steps = grad.getSteps();
        if (steps > 0) {
            order = grad.getOrder();
            amplitudeArray = new double[steps];
            for (int i = 0; i < steps; i++) {
                amplitudeArray[i] = -grad.getAmplitudeArray(i) * grad.getEquivalentTime() / equivalentTime;
            }
        }
        return amplitudeArray;
    }

    public void reoderPhaseEncoding(TransformPlugin plugin, int echoTrainLength, int acquisitionMatrixDimension2D, int acquisitionMatrixDimension1D) {
        double loopNumber, indexNew;
        if (amplitudeArray != null) {
            double[] newTable = new double[acquisitionMatrixDimension2D];
            for (int j = 0; j < acquisitionMatrixDimension2D; j++) {
                int[] indexScan = plugin.invTransf(0, j, 0, 0);
                loopNumber = indexScan[0] / acquisitionMatrixDimension1D; // Echo-block number: ETL-loop index
                indexNew = indexScan[1] * echoTrainLength + loopNumber;    // indexScan[1]: index de Nb 2D
                newTable[(int) indexNew] = amplitudeArray[j];
            }
            amplitudeArray = newTable;
        }
    }


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  Spoiler
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    public boolean addSpoiler(double gradAmplitude) {
        boolean testSpoilerSupThan100 = true;
        if (Double.isNaN(amplitude)) {
            amplitude = gradAmplitude;
        } else {
            amplitude += gradAmplitude;
        }
        calculateStaticArea();
        double[] gradMaxMin = checkGradientMax();
        if (gradMaxMin[0] > 100.0) {
            amplitude = 100.0;
            spoilerExcess = gradMaxMin[0] - 100.0;
            testSpoilerSupThan100 = false;
        }
        return (testSpoilerSupThan100);
    }


    // ---------------------------------------------------------------
    // ----------------- General Methode----------------------------------------------
    private double getInferiorSpectralWidth(double spectral_width) {
        return 3906250.0 / (Math.ceil(3906250.0 / spectral_width) + 1);
    }

    private void setSequenceTableValues(Table table, Order order, double... values) {
        table.clear();
        table.setOrder(order);
        table.setLocked(true);
        for (double value : values) {
            table.add(value);
        }
    }

}