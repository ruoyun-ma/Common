package common;

import rs2d.spinlab.instrument.util.GradientMath;
import rs2d.spinlab.sequence.Sequence;
import rs2d.spinlab.sequence.table.Shape;
import rs2d.spinlab.sequence.table.Table;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.utility.Nucleus;

/**
 * Commons V2.8: add instrument interface
 * Class Gradient
 * V2.9 with ratio for non proton
 * V2.3- 2019-06-06 JR from DW EPI
 * V2.1- 2018-03-20b JR
 */
public class Gradient2Event extends Gradient {

    protected Gradient2Event gradFlowComp = null;

    protected static double gMax = GradientMath.getMaxGradientStrength();

    public Gradient2Event(Table amplitudeTab, Shape shapeUpTab, Shape shapeDownTab, Table rampTimeUpTab, Table rampTimeDownTab) {
        super(amplitudeTab, null, shapeUpTab, shapeDownTab, rampTimeUpTab, rampTimeDownTab);
    }


    public static Gradient2Event createGradient(Sequence sequence, GeneratorSequenceParamEnum amplitudeTab, GeneratorSequenceParamEnum shapeUpTab, GeneratorSequenceParamEnum shapeDownTab, GeneratorSequenceParamEnum rampTimeTab) {
        return new Gradient2Event(sequence.getPublicTable(amplitudeTab.name()), (Shape) sequence.getPublicTable(shapeUpTab.name()),
                (Shape) sequence.getPublicTable(shapeDownTab.name()), sequence.getPublicTable(rampTimeTab.name()), sequence.getPublicTable(rampTimeTab.name()));
    }

    public static Gradient2Event createGradient(Sequence sequence, GeneratorSequenceParamEnum amplitudeTab, GeneratorSequenceParamEnum shapeUpTab, GeneratorSequenceParamEnum shapeDownTab, GeneratorSequenceParamEnum rampTimeUpTab, GeneratorSequenceParamEnum rampTimeDownTab) {
        return new Gradient2Event(sequence.getPublicTable(amplitudeTab.name()), (Shape) sequence.getPublicTable(shapeUpTab.name()),
                (Shape) sequence.getPublicTable(shapeDownTab.name()), sequence.getPublicTable(rampTimeUpTab.name()), sequence.getPublicTable(rampTimeDownTab.name()));
    }


//    For non proton
    public Gradient2Event(Table amplitudeTab, Shape shapeUpTab, Shape shapeDownTab, Table rampTimeUpTab, Table rampTimeDownTab, Nucleus nucleus) {
        super(amplitudeTab, null, shapeUpTab, shapeDownTab, rampTimeUpTab, rampTimeDownTab,nucleus);
    }



    public static Gradient2Event createGradient(Sequence sequence, GeneratorSequenceParamEnum amplitudeTab, GeneratorSequenceParamEnum shapeUpTab, GeneratorSequenceParamEnum shapeDownTab, GeneratorSequenceParamEnum rampTimeTab, Nucleus nucleus)  {
        return new Gradient2Event(sequence.getPublicTable(amplitudeTab.name()), (Shape) sequence.getPublicTable(shapeUpTab.name()),
                (Shape) sequence.getPublicTable(shapeDownTab.name()), sequence.getPublicTable(rampTimeTab.name()), sequence.getPublicTable(rampTimeTab.name()),nucleus);
    }

    public static Gradient2Event createGradient(Sequence sequence, GeneratorSequenceParamEnum amplitudeTab, GeneratorSequenceParamEnum shapeUpTab, GeneratorSequenceParamEnum shapeDownTab, GeneratorSequenceParamEnum rampTimeUpTab, GeneratorSequenceParamEnum rampTimeDownTab, Nucleus nucleus)  {
        return new Gradient2Event(sequence.getPublicTable(amplitudeTab.name()), (Shape) sequence.getPublicTable(shapeUpTab.name()),
                (Shape) sequence.getPublicTable(shapeDownTab.name()), sequence.getPublicTable(rampTimeUpTab.name()), sequence.getPublicTable(rampTimeDownTab.name()),nucleus);
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  general  methodes
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     * Calculate equivalentTime of a rectangular gradient with same Area and amplitude
     *
     */
    @Override
    public void prepareEquivalentTime() {
        if (Double.isNaN(grad_shape_rise_time)) {
            computeShapeRiseTime();
        }
        equivalentTime = grad_shape_rise_time;
    }

    @Override
    public void refocalizeGradient() {
        calculateStaticAmplitude();
    }

    public void refocalizeGradient(Gradient2Event grad, double ratio) {
        bStaticGradient = true;
        staticArea = -grad.getStaticArea() * ratio;
        calculateStaticAmplitude();
    }

    public void refocalizeGradientWithFlowComp(Gradient2Event grad, double ratio, Gradient2Event gradflowcomp) {
        gradFlowComp = gradflowcomp;
        // to modify , flow Comp
        //to do: modify the calculation and prepare as well gradFlowComp Gradient


        bStaticGradient = true;
        staticArea = -grad.getStaticArea() * ratio;
        calculateStaticAmplitude();
    }

    public boolean refocalizeGradientWithAmplitude(Gradient2Event grad, double ratio, double amplitude) {
        if (Double.isNaN(grad_shape_rise_time)) {
            computeShapeRiseTime();
        }
        staticArea = -grad.getStaticArea() * ratio;
        boolean test_Amplitude = true;
        equivalentTime = staticArea / amplitude;
        double topTime = equivalentTime - grad_shape_rise_time;
        if (topTime < 0.000004) {
            topTime = 0.000004;
//            flatTimeTable.set(0, topTime);
//            prepareEquivalentTime();
//            calculateStaticAmplitude();
            test_Amplitude = false;
        }
        prepareEquivalentTime();
        calculateStaticAmplitude();
        return test_Amplitude;
    }


    /**
     * calculate READOUT refocusing gradient Amplitude handeling ETL
     *
     * @param grad : Readout Gradient
     */
    public void refocalizeReadoutGradient(Gradient2Event grad, double ratio) {
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
    public void refocalizeReadoutGradients(Gradient2Event grad, double ratio) {
        steps = grad.getSteps();
        amplitudeArray = new double[steps];
        for (int i = 0; i < steps; i++) {
            // flatTimeTable.get(0).doubleValue() + grad_shape_rise_time
            amplitudeArray[i] = -grad.getAmplitudeArray(i) * grad.getEquivalentTime() / this.getEquivalentTime() * ratio;
        }

        order = grad.getOrder();
    }


    public void preparePhaseEncodingForCheckWithFlowComp(int matrixDimensionForCheck, int matrixDimension, double fovDim, boolean isKSCentred, Gradient2Event gradflowcomp) {
        gradFlowComp = gradflowcomp;
        // to modify , flow Comp
        //to do: modify the calculation and prepare as well gradFlowComp Gradient

        double grad_total_area_phase = prepPhaseGradTotalArea(matrixDimensionForCheck, fovDim);
        double grad_index_max_phase = prepPhaseGradIndexMax(isKSCentred);
        maxAreaPE = grad_index_max_phase * grad_total_area_phase;

        preparePhaseEncoding(matrixDimension, fovDim, isKSCentred);
    }

    public double[] refocalizePhaseEncodingGradient(Gradient2Event grad) {
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

}