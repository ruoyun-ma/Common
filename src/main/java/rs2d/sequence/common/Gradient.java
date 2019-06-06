package rs2d.sequence.common;

import rs2d.spinlab.data.transformPlugin.TransformPlugin;
import rs2d.spinlab.hardware.controller.HardwareHandler;
import rs2d.spinlab.hardware.devices.DeviceManager;
import rs2d.spinlab.instrument.util.GradientMath;
import rs2d.spinlab.sequence.Sequence;
import rs2d.spinlab.sequence.table.Shape;
import rs2d.spinlab.sequence.table.Table;
import rs2d.spinlab.sequence.table.Utility;
import rs2d.spinlab.tools.table.Order;

/**
 * Class Gradient
 * V2.3- 2019-06-06 JR from DW EPI
 * V2.2- 2018-12-19 JR
 * V2.1- 2017-10-24 JR
 */
public class Gradient {
    protected Table amplitudeTable = null;
    protected Table flatTimeTable = null;
    protected Shape shapeUpTable = null;
    protected Shape shapeDownTable = null;
    protected Table rampTimeUpTable;
    protected Table rampTimeDownTable;

    protected double amplitude = Double.NaN;
    protected double staticArea = Double.NaN;

    protected double[] amplitudeArray = null;
    protected double maxAreaPE = Double.NaN;

    protected int steps = -1;
    protected Order order = Order.FourLoop;

    protected double grad_shape_rise_time = Double.NaN;
    protected double equivalentTime = Double.NaN;

    // Slc
    protected double sliceThicknessExcitation = Double.NaN;
    protected double txBandwidth = Double.NaN;

    // RO
    protected double spectralWidth = Double.NaN;

    // PE
    protected boolean bPhaseEncoding = false;
    protected double fovPhase = Double.NaN;
    protected boolean isKSCentred = false;

    protected double k0pos = Double.NaN;
    protected double spoilerExcess = Double.NaN;
    protected double minTopTime = Double.NaN;

    protected boolean bRefocalizeGradient = false;
    protected boolean bStaticGradient = false;

    protected Gradient gradFlowComp = null;

    protected static double gMax = GradientMath.getMaxGradientStrength();

    public Gradient(Table amplitudeTab, Table flat_TimeTab, Shape shapeUpTab, Shape shapeDownTab, Table rampTimeUpTab, Table rampTimeDownTab) {
        amplitudeTable = amplitudeTab;
        flatTimeTable = flat_TimeTab;
        shapeUpTable = shapeUpTab;
        shapeDownTable = shapeDownTab;
        rampTimeUpTable = rampTimeUpTab;
        rampTimeDownTable = rampTimeDownTab;
        gMax = GradientMath.getMaxGradientStrength();
        init();
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

    protected void init() {
        computeShapeRiseTime();
        prepareEquivalentTime();
    }

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

    public double getAmplitudeArray_mTpm(int pos) {
        if (amplitudeArray == null || pos >= steps) {
            return Double.NaN;
        } else {
            return amplitudeArray[pos] * gMax / 100.0;
        }
    }

    public double getSliceThickness() {
        return sliceThicknessExcitation;
    }

    public double getSpoilerExcess() {
        return spoilerExcess;
    }

    public double getMinTopTime() {
        return minTopTime;
    }

    public int getSteps() {
        return steps;
    }

    public Order getOrder() {
        return order;
    }

    public double getK0pos() {
        return k0pos;
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

    public double getTotalAbsArea() {
        return Math.abs(getTotalArea());
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

    public void setAmplitude(double... values) {
        if (values.length == 1) {
            amplitude = values[0];
        } else {
            amplitudeArray = new double[values.length];
            int i = 0;
            for (double value : values) {
                //System.out.println(i + " " + value);
                amplitudeArray[i] = value;
                i += 1;
            }
            steps = i;
        }
    }

    public void setAmplitude(Order order, double... values) {
        setAmplitude(values);
        applyAmplitude(order);
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
    protected double computeShapeRiseTime() {
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
        minTopTime = flatTimeTable.get(0).doubleValue();
        equivalentTime = (flatTimeTable.get(0).doubleValue() + grad_shape_rise_time);
        return equivalentTime;
    }

    public void refocalizeGradient() {
        calculateStaticAmplitude();
    }

    public void refocalizeGradient(Gradient grad, double ratio) {
        bStaticGradient = true;
        staticArea = -grad.getStaticArea() * ratio;
        calculateStaticAmplitude();
    }

    public void refocalizeGradientWithFlowComp(Gradient grad, double ratio, Gradient gradflowcomp) {
        gradFlowComp = gradflowcomp;
        // to modify , flow Comp
        //to do: modify the calculation and prepare as well gradFlowComp Gradient


        bStaticGradient = true;
        staticArea = -grad.getStaticArea() * ratio;
        calculateStaticAmplitude();
    }

    public boolean refocalizeGradientWithAmplitude(Gradient grad, double ratio, double amplitude) {
        if (grad_shape_rise_time == Double.NaN) {
            computeShapeRiseTime();
        }
        staticArea = -grad.getStaticArea() * ratio;
        boolean test_Amplitude = true;
        equivalentTime = staticArea / amplitude;
        double topTime = equivalentTime - grad_shape_rise_time;
        if (topTime < 0.000004) {
            topTime = 0.000004;
            test_Amplitude = false;
        }
        flatTimeTable.set(0, topTime);
        prepareEquivalentTime();
        calculateStaticAmplitude();
        return test_Amplitude;
    }

    public void rePrepare() {
        prepareEquivalentTime();
        if (bPhaseEncoding)
            preparePhaseEncoding();
        if (bStaticGradient)
            calculateStaticAmplitude();

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
        double spectralWidth_init = (gMax * GradientMath.GAMMA * fov);
        double spectralWidth = spectralWidth_init;
        while (DeviceManager.getInstance().getCompilerLoader().getCompiler().getSequenceCompiler().getNearestSW(spectralWidth) >= spectralWidth_init) {
            // Hardware.getSequenceCompiler() //future version
            spectralWidth = getInferiorSpectralWidth(spectralWidth);
        }
        spectralWidth = DeviceManager.getInstance().getCompilerLoader().getCompiler().getSequenceCompiler().getNearestSW(spectralWidth);
        // Hardware.getSequenceCompiler() //future version
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

    /*
     * calculate READOUT refocusing
     *
     * @param grad : Readout Gradient
     * @param ratio : ratio to compensate
     */
    public void refocalizeReadoutGradients(Gradient grad, double ratio) {
        steps = grad.getSteps();
        amplitudeArray = new double[steps];
        for (int i = 0; i < steps; i++) {
            // flatTimeTable.get(0).doubleValue() + grad_shape_rise_time
            amplitudeArray[i] = -grad.getAmplitudeArray(i) * grad.getEquivalentTime() / this.getEquivalentTime() * ratio;
        }

        order = grad.getOrder();
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
            sliceThicknessExcitation = ceilToSubDecimal(tx_bandwidth / ((GradientMath.GAMMA) * gMax), 6);
            amplitude = (tx_bandwidth / ((GradientMath.GAMMA) * sliceThicknessExcitation)) * 100.0 / gMax;                 // amplitude in T/m
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

    public boolean prepareEPIBlip(int step, double fovDim) {
        staticArea = prepPhaseGradTotalArea(step, fovDim);
        amplitude = staticArea / equivalentTime;
        return amplitude <= 100.0;
    }

    public double getShapeFactor() {
        double grad_shape_rise_factor_up = Utility.voltageFillingFactor(shapeUpTable);
        double grad_shape_rise_factor_down = Utility.voltageFillingFactor(shapeDownTable);
        return (grad_shape_rise_factor_up + grad_shape_rise_factor_down) / 2;
    }

    public void preparePhaseEncodingForCheck(int matrixDimensionForCheck, int matrixDimension, double fovDim, boolean isKSCentred) {
        double grad_total_area_phase = prepPhaseGradTotalArea(matrixDimensionForCheck, fovDim);
        double grad_index_max_phase = prepPhaseGradIndexMax(isKSCentred);
        maxAreaPE = grad_index_max_phase * grad_total_area_phase;
        preparePhaseEncoding(matrixDimension, fovDim, isKSCentred);
    }


    public void preparePhaseEncodingForCheckWithFlowComp(int matrixDimensionForCheck, int matrixDimension, double fovDim, boolean isKSCentred, Gradient gradflowcomp) {
        gradFlowComp = gradflowcomp;
        // to modify , flow Comp
        //to do: modify the calculation and prepare as well gradFlowComp Gradient

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
        k0pos = gradIndexMaxPhase * steps;
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

    public double[] refocalizePhaseEncodingGradientSEEPI(Gradient grad, Gradient gradBlip, int echoTrainLength) {
        steps = grad.getSteps();
        System.out.println("      ");
        System.out.println("      ");
        System.out.println("      ");
        System.out.println("steps      " + steps);

        if (steps > 0) {
            order = grad.getOrder();
            amplitudeArray = new double[steps];
            for (int i = 0; i < steps; i++) {
                amplitudeArray[i] = -(-grad.getAmplitudeArray(i) * grad.getEquivalentTime()
                        + gradBlip.getAmplitude() * gradBlip.getEquivalentTime() * echoTrainLength) / equivalentTime;
            }
            System.out.println("gradBlip.getAmplitude()  " + gradBlip.getAmplitude());
            System.out.println("gradBlip.getEquivalentTime()  " + gradBlip.getEquivalentTime());
        }
        return amplitudeArray;
    }


    public void reoderPhaseEncodingForSEEPI(int echoTrainLength) {
        // flow Comp
        int new_steps = steps / echoTrainLength;
        double[] newTable = new double[new_steps];
        int fact = 1;
        for (int i = 0; i < new_steps; i++) {
            int ii = i * fact;
            newTable[i] = -amplitudeArray[ii];
        }
        amplitudeArray = newTable;
        steps = new_steps;
    }


    public void reoderPhaseEncoding(TransformPlugin plugin, int echoTrainLength, int acquisitionMatrixDimension2D, int acquisitionMatrixDimension1D) {
        // flow Comp
        if (gradFlowComp != null) {
            gradFlowComp.reoderPhaseEncoding(plugin, echoTrainLength, acquisitionMatrixDimension2D, acquisitionMatrixDimension1D);
        }
        double loopNumber, indexNew;
        if (amplitudeArray != null) {
            double[] newTable = new double[acquisitionMatrixDimension2D];
            int[] tmp = Centric(acquisitionMatrixDimension2D);
            for (int j = 0; j < acquisitionMatrixDimension2D; j++) {
                int[] indexScan = plugin.invTransf(0, j, 0, 0);
                 if ("Centric4D".equalsIgnoreCase(plugin.getName())) {
                    indexScan[1] = tmp[j];
                }
                loopNumber = indexScan[0] / acquisitionMatrixDimension1D; // Echo-block number: ETL-loop index
                indexNew = indexScan[1] * echoTrainLength + loopNumber;    // indexScan[1]: index de Nb 2D
                newTable[(int) indexNew] = amplitudeArray[j];
            }
            amplitudeArray = newTable;
        }
    }

    public void reoderPhaseEncoding3D(TransformPlugin plugin, int acquisitionMatrixDimension3D) {
        // to modify , flow Comp
        if (gradFlowComp != null) {
            gradFlowComp.reoderPhaseEncoding3D(plugin, acquisitionMatrixDimension3D);
        }

        double indexNew;
        if (amplitudeArray != null) {
            double[] newTable = new double[acquisitionMatrixDimension3D];

            System.out.println("----- " + plugin.getName());
            int[] tmp = Centric(acquisitionMatrixDimension3D);
            for (int k = 0; k < acquisitionMatrixDimension3D; k++) {
                int[] indexScan = plugin.invTransf(0, 0, k, 0);
                //              System.out.println(j + " :  " +indexScan[0] + " " +indexScan[1]+ "        " +  plugin.transf(0, j, 0, 0)[0]+ " " +  plugin.transf(0, j, 0, 0)[1]);
                if ("Centric4D".equalsIgnoreCase(plugin.getName())) {
                    indexScan[2] = tmp[k];
//                    System.out.println(j + " :  " + indexScan[1]);
                }
                indexNew = indexScan[2];    // indexScan[1]: index de Nb 2D
                newTable[(int) indexNew] = amplitudeArray[k];
            }
            amplitudeArray = newTable;
        }
    }


    public int[] Centric(int acquisition_matrix_dimension_3D) {
        int[] k_index = new int[acquisition_matrix_dimension_3D];
        int[] tmpInv = new int[acquisition_matrix_dimension_3D];

        k_index[0] = 0;
        for (int i = 1; i < acquisition_matrix_dimension_3D; i++) {
            if (k_index[i - 1] == 0) {
                k_index[i] = 1;
            } else if (k_index[i - 1] > 0) {
                k_index[i] = -1 * k_index[i - 1];
            } else if (k_index[i - 1] < 0) {
                k_index[i] = -1 * (k_index[i - 1] - 1);
            }
        }
        for (int i = 0; i < acquisition_matrix_dimension_3D; i++) {
            k_index[i] = k_index[i] + (acquisition_matrix_dimension_3D / 2 - 1);
        }
        for (int i = 0; i < acquisition_matrix_dimension_3D; i++) {
            tmpInv[k_index[i]] = i;
        }
        return tmpInv;
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  Spoiler
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    public boolean addSpoiler(double gradAmplitude) {
        boolean testSpoilerSupThan100 = true;
        bStaticGradient = true;
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
            minTopTime = ceilToSubDecimal((gradMaxMin[0] * equivalentTime - grad_shape_rise_time * 100.0) / 100.0, 5);
            testSpoilerSupThan100 = false;
        }
        return (testSpoilerSupThan100);
    }

    public boolean addSpoiler(double pixel_dimension, double factor) {
        bStaticGradient = true;
        double grad_area_spoiler = factor / ((GradientMath.GAMMA) * pixel_dimension);//GradientMath.GAMMA: gamma/2pi értéke Hz/T-ban
        double grad_amp_spoiler = (grad_area_spoiler / equivalentTime) / gMax * 100.0;//
        return (addSpoiler(grad_amp_spoiler));
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  Flow Compensation
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    public double getMoment(double t, String type, int order) {
        double rise_time_up = rampTimeUpTable.get(0).doubleValue();
        double rise_time_down = rampTimeDownTable.get(0).doubleValue();
        double plato = flatTimeTable.get(0).doubleValue();
        double amp = amplitude;
        if (order == 0) {

            if (type == "slice") {
                t = 0;
                double moment0 = amp * (2 * (rise_time_up + rise_time_down) / Math.PI + plato) / 2;
                return moment0;
            } else if (type == "read") {
                double moment0 = amp * (2 * (rise_time_up + rise_time_down) / Math.PI + plato) / 2;
                return moment0;
            } else {
                double moment0 = amp * (2 * (rise_time_up + rise_time_down) / Math.PI + plato);
                return moment0;
            }

        } else if (type == "slice") {
            t = 0;
            double moment0 = amp * (2 * (rise_time_up + rise_time_down) / Math.PI + plato) / 2;
            double moment1 = 1 / 8 * amp * plato * plato + amp * (rise_time_up + rise_time_down) / 2 * (Math.PI * plato + (Math.PI - 2.0) * (rise_time_up + rise_time_down)) / Math.PI / Math.PI;
            double moment1_all = moment1 + t * moment0;
            return moment1_all;
        } else if (type == "read") {
            double moment0 = amp * (2 * (rise_time_up + rise_time_down) / Math.PI + plato) / 2;
            double moment1 = 1 / 2 * amp * ((rise_time_up + rise_time_down) / 2 * plato + plato * plato) + amp * (rise_time_up + rise_time_down) * (rise_time_up + rise_time_down) / Math.PI / Math.PI;
            double moment1_all = moment1 + t * moment0;
            return moment1_all;
        } else {
            double moment0 = amp * (2 * (rise_time_up + rise_time_down) / Math.PI + plato);
            double moment1 = amp * (plato + (rise_time_up + rise_time_down)) * ((rise_time_up + rise_time_down) / Math.PI + plato / 2);
            double moment1_all = moment1 + t * moment0;
            return moment1_all;
        }
    }

    // ---------------------------------------------------------------
    // ----------------- General Methode----------------------------------------------
    private double getInferiorSpectralWidth(double spectral_width)throws Exception {
        double  SWmin = DeviceManager.getInstance().getCompilerLoader().getCompiler().getSequenceCompiler().getMaxSW();
        double SWmax = Math.round(DeviceManager.getInstance().getCompilerLoader().getCompiler().getSequenceCompiler().getMaxSW());
        // Hardware.getSequenceCompiler() //future version
//        return 3906250.0 / (Math.ceil(3906250.0 / spectral_width) + 1);
        return SWmax / (Math.ceil(SWmax / spectral_width) + 1);
    }

    private double ceilToSubDecimal(double numberToBeRounded, double Order) {
        return Math.ceil(numberToBeRounded * Math.pow(10, Order)) / Math.pow(10, Order);
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