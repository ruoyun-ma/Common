package rs2d.sequence.common;

import rs2d.spinlab.data.transformPlugin.TransformPlugin;
import rs2d.spinlab.hardware.controller.HardwareHandler;
import rs2d.spinlab.instrument.util.GradientMath;
import rs2d.spinlab.sequence.Sequence;
import rs2d.spinlab.sequence.table.Shape;
import rs2d.spinlab.sequence.table.Table;
import rs2d.spinlab.sequence.table.Utility;
import rs2d.spinlab.sequenceGenerator.GeneratorSequenceParamEnum;
import rs2d.spinlab.tools.table.Order;

/**
 * Class Gradient
 * V2.3- 2019-06-06 JR from DW EPI
 * V2.1- 2018-03-20b JR
 */
public class Gradient5Event extends Gradient {

    protected Table flatTimeTable1 = null;
    protected Table flatTimeTable3 = null;

    public Gradient5Event(Table amplitudeTab, Table flat_TimeTab1, Table flat_TimeTab2, Table flat_TimeTab3, Shape shapeUpTab, Shape shapeDownTab, Table rampTimeUpTab, Table rampTimeDownTab) {
        super(amplitudeTab, flat_TimeTab2, shapeUpTab, shapeDownTab, rampTimeUpTab, rampTimeDownTab);
        flatTimeTable1 = flat_TimeTab1;
        flatTimeTable3 = flat_TimeTab3;
        init();
        //        Table flatTimeTable1 = getSequence().getPublicTable(Time_grad_read_crusher);
//        System.out.println(flatTimeTable1.getFirst());
//        System.out.println(flatTimeTable1.size());
//        System.out.println(flatTimeTable1.getName());
//        System.out.println(flatTimeTable1.get(0).doubleValue());

    }

    public static Gradient5Event createGradient(Sequence sequence, GeneratorSequenceParamEnum amplitudeTab, GeneratorSequenceParamEnum flat_TimeTab1, GeneratorSequenceParamEnum flat_TimeTab2, GeneratorSequenceParamEnum flat_TimeTab3, GeneratorSequenceParamEnum shapeUpTab, GeneratorSequenceParamEnum shapeDownTab, GeneratorSequenceParamEnum rampTimeTab) {
//        System.out.println(" ");
//        System.out.println(" Gradient5Event createGradient 1  ");
//        System.out.println(" flat_TimeTab1  " + flat_TimeTab1);
//        System.out.println(flat_TimeTab1);
        return new Gradient5Event(sequence.getPublicTable(amplitudeTab.name()), sequence.getPublicTable(flat_TimeTab1.name()), sequence.getPublicTable(flat_TimeTab2.name()), sequence.getPublicTable(flat_TimeTab3.name()), (Shape) sequence.getPublicTable(shapeUpTab.name()),
                (Shape) sequence.getPublicTable(shapeDownTab.name()), sequence.getPublicTable(rampTimeTab.name()), sequence.getPublicTable(rampTimeTab.name()));
    }

    public static Gradient5Event createGradient(Sequence sequence, GeneratorSequenceParamEnum amplitudeTab, GeneratorSequenceParamEnum flat_TimeTab1, GeneratorSequenceParamEnum flat_TimeTab2, GeneratorSequenceParamEnum flat_TimeTab3, GeneratorSequenceParamEnum shapeUpTab, GeneratorSequenceParamEnum shapeDownTab, GeneratorSequenceParamEnum rampTimeUpTab, GeneratorSequenceParamEnum rampTimeDownTab) {
//        System.out.println(" ");
//        System.out.println(" Gradient5Event createGradient 2 ");
//        System.out.println(" flat_TimeTab1  " + flat_TimeTab1);
        return new Gradient5Event(sequence.getPublicTable(amplitudeTab.name()), sequence.getPublicTable(flat_TimeTab1.name()), sequence.getPublicTable(flat_TimeTab2.name()), sequence.getPublicTable(flat_TimeTab3.name()), (Shape) sequence.getPublicTable(shapeUpTab.name()),
                (Shape) sequence.getPublicTable(shapeDownTab.name()), sequence.getPublicTable(rampTimeUpTab.name()), sequence.getPublicTable(rampTimeDownTab.name()));
    }


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //                  general  methodes
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    protected void init() {
        computeShapeRiseTime();
        prepareEquivalentTime();
    }


    /**
     * Calculate equivalentTime of a rectangular gradient with same Area and amplitude
     *
     * @return equivalentTime :
     *//*
     */
    @Override
    public double prepareEquivalentTime() {
//        System.out.println(" ");
//        System.out.println(" prepareEquivalentTime() ");
//        System.out.println(" flatTimeTable " + flatTimeTable1);

        if (grad_shape_rise_time == Double.NaN) {
            computeShapeRiseTime();
        }
        if (flatTimeTable1  != null) {
//            System.out.println("flatTimeTable.get(0).doubleValue()  = " + flatTimeTable.get(0).doubleValue());
//            System.out.println("flatTimeTable1.get(0).doubleValue()  = " + flatTimeTable1.get(0).doubleValue());
//            System.out.println("flatTimeTable3.get(0).doubleValue()  = " + flatTimeTable3.get(0).doubleValue());
            minTopTime = flatTimeTable.get(0).doubleValue() + flatTimeTable1.get(0).doubleValue() + flatTimeTable3.get(0).doubleValue();
            equivalentTime = (minTopTime + grad_shape_rise_time);
//            System.out.println("minTopTime  = " + minTopTime);
//            System.out.println("grad_shape_rise_time  = " + grad_shape_rise_time);
        }
        return equivalentTime;
    }
/*

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //     RO             Readout  methodes             RO
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    */

    /**
     * calculate readout Gradient Amplitude
     *
     * @param spectralWidth
     * @param fov
     * @return testSpectralWidth : false if the Spectralwidth need to be nicreased(call getSpectralWidth() )
     *//*
     */
    @Override
    public boolean calculateReadoutGradient(double spectralWidth, double fov) throws Exception {
        boolean testSpectralWidth = super.calculateReadoutGradient(spectralWidth, fov);
        calculateStaticArea();
        return testSpectralWidth;
    }


}