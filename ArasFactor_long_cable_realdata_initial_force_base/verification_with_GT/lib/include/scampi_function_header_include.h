// ************************************ model_include/FK ************************************
// cost1 ---------
#include "model_include/FK_residual_func_cost1.h"
#include "model_include/FK_residual_func_cost1_wrt_fh1.h"
#include "model_include/FK_residual_func_cost1_wrt_fv1.h"
#include "model_include/FK_residual_func_cost1_wrt_position_vector.h"
#include "model_include/FK_residual_func_cost1_wrt_DeltaRot.h"

// cost2 ---------
#include "model_include/FK_residual_func_cost2.h"
#include "model_include/FK_residual_func_cost2_wrt_fh1.h"
#include "model_include/FK_residual_func_cost2_wrt_fv1.h"
#include "model_include/FK_residual_func_cost2_wrt_position_vector.h"
#include "model_include/FK_residual_func_cost2_wrt_DeltaRot.h"
// cost3 ---------
#include "model_include/FK_residual_func_cost3.h"
#include "model_include/FK_residual_func_cost3_wrt_fh1.h"
#include "model_include/FK_residual_func_cost3_wrt_fv1.h"
#include "model_include/FK_residual_func_cost3_wrt_position_vector.h"
#include "model_include/FK_residual_func_cost3_wrt_DeltaRot.h"


// ************************************ IK ************************************
// cost1 ---------
// NL0 ------------
#include "model_include/IK_residual_func_cost1_Nl0.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl0.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl0.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl0.h"
// cost2 ---------
// Nl0 ------------
#include "model_include/IK_residual_func_cost2_Nl0.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl0.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl0.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl0.h"
// cost3 ---------
// Nl0 ------------
#include "model_include/IK_residual_func_cost3_Nl0.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl0.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl0.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl0.h"

// cost1 ---------
// NL1 ------------
#include "model_include/IK_residual_func_cost1_Nl1.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl1.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl1.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl1.h"
// cost2 ---------
// Nl1 ------------
#include "model_include/IK_residual_func_cost2_Nl1.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl1.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl1.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl1.h"
// cost3 ---------
// Nl1 ------------
#include "model_include/IK_residual_func_cost3_Nl1.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl1.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl1.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl1.h"

// cost1 ---------
// NL2 ------------
#include "model_include/IK_residual_func_cost1_Nl2.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl2.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl2.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl2.h"
// cost2 ---------
// Nl2 ------------
#include "model_include/IK_residual_func_cost2_Nl2.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl2.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl2.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl2.h"
// cost3 ---------
// Nl2 ------------
#include "model_include/IK_residual_func_cost3_Nl2.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl2.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl2.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl2.h"

// cost1 ---------
// NL3 ------------
#include "model_include/IK_residual_func_cost1_Nl3.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl3.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl3.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl3.h"
// cost2 ---------
// Nl3 ------------
#include "model_include/IK_residual_func_cost2_Nl3.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl3.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl3.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl3.h"
// cost3 ---------
// Nl3 ------------
#include "model_include/IK_residual_func_cost3_Nl3.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl3.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl3.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl3.h"

// cost1 ---------
// NL4 ------------
#include "model_include/IK_residual_func_cost1_Nl4.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl4.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl4.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl4.h"
// cost2 ---------
// Nl4 ------------
#include "model_include/IK_residual_func_cost2_Nl4.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl4.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl4.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl4.h"
// cost3 ---------
// Nl4 ------------
#include "model_include/IK_residual_func_cost3_Nl4.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl4.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl4.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl4.h"

// cost1 ---------
// NL5 ------------
#include "model_include/IK_residual_func_cost1_Nl5.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl5.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl5.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl5.h"
// cost2 ---------
// Nl5 ------------
#include "model_include/IK_residual_func_cost2_Nl5.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl5.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl5.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl5.h"
// cost3 ---------
// Nl5 ------------
#include "model_include/IK_residual_func_cost3_Nl5.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl5.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl5.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl5.h"

// cost1 ---------
// NL6 ------------
#include "model_include/IK_residual_func_cost1_Nl6.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl6.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl6.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl6.h"
// cost2 ---------
// Nl6 ------------
#include "model_include/IK_residual_func_cost2_Nl6.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl6.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl6.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl6.h"
// cost3 ---------
// Nl6 ------------
#include "model_include/IK_residual_func_cost3_Nl6.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl6.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl6.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl6.h"

// cost1 ---------
// NL7 ------------
#include "model_include/IK_residual_func_cost1_Nl7.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl7.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl7.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl7.h"
// cost2 ---------
// Nl7 ------------
#include "model_include/IK_residual_func_cost2_Nl7.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl7.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl7.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl7.h"
// cost3 ---------
// Nl7 ------------
#include "model_include/IK_residual_func_cost3_Nl7.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl7.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl7.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl7.h"

// cost1 ---------
// NL8 ------------
#include "model_include/IK_residual_func_cost1_Nl8.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl8.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl8.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl8.h"
// cost2 ---------
// Nl8 ------------
#include "model_include/IK_residual_func_cost2_Nl8.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl8.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl8.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl8.h"
// cost3 ---------
// Nl8 ------------
#include "model_include/IK_residual_func_cost3_Nl8.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl8.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl8.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl8.h"

// cost1 ---------
// NL9 ------------
#include "model_include/IK_residual_func_cost1_Nl9.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl9.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl9.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl9.h"
// cost2 ---------
// Nl9 ------------
#include "model_include/IK_residual_func_cost2_Nl9.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl9.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl9.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl9.h"
// cost3 ---------
// Nl9 ------------
#include "model_include/IK_residual_func_cost3_Nl9.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl9.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl9.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl9.h"

// cost1 ---------
// NL10 ------------
#include "model_include/IK_residual_func_cost1_Nl10.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl10.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl10.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl10.h"


// cost2 ---------
// Nl10 ------------
#include "model_include/IK_residual_func_cost2_Nl10.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl10.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl10.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl10.h"


// cost3 ---------
// Nl10 ------------
#include "model_include/IK_residual_func_cost3_Nl10.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl10.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl10.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl10.h"



// cost1 ---------
// NL11 ------------
#include "model_include/IK_residual_func_cost1_Nl11.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl11.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl11.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl11.h"


// cost2 ---------
// Nl11 ------------
#include "model_include/IK_residual_func_cost2_Nl11.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl11.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl11.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl11.h"


// cost3 ---------
// Nl11 ------------
#include "model_include/IK_residual_func_cost3_Nl11.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl11.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl11.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl11.h"



// cost1 ---------
// NL12 ------------
#include "model_include/IK_residual_func_cost1_Nl12.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl12.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl12.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl12.h"


// cost2 ---------
// Nl12 ------------
#include "model_include/IK_residual_func_cost2_Nl12.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl12.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl12.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl12.h"


// cost3 ---------
// Nl12 ------------
#include "model_include/IK_residual_func_cost3_Nl12.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl12.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl12.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl12.h"



// cost1 ---------
// NL13 ------------
#include "model_include/IK_residual_func_cost1_Nl13.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl13.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl13.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl13.h"


// cost2 ---------
// Nl13 ------------
#include "model_include/IK_residual_func_cost2_Nl13.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl13.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl13.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl13.h"


// cost3 ---------
// Nl13 ------------
#include "model_include/IK_residual_func_cost3_Nl13.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl13.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl13.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl13.h"



// cost1 ---------
// NL14 ------------
#include "model_include/IK_residual_func_cost1_Nl14.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl14.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl14.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl14.h"


// cost2 ---------
// Nl14 ------------
#include "model_include/IK_residual_func_cost2_Nl14.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl14.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl14.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl14.h"


// cost3 ---------
// Nl14 ------------
#include "model_include/IK_residual_func_cost3_Nl14.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl14.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl14.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl14.h"



// cost1 ---------
// NL15 ------------
#include "model_include/IK_residual_func_cost1_Nl15.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl15.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl15.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl15.h"


// cost2 ---------
// Nl15 ------------
#include "model_include/IK_residual_func_cost2_Nl15.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl15.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl15.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl15.h"


// cost3 ---------
// Nl15 ------------
#include "model_include/IK_residual_func_cost3_Nl15.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl15.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl15.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl15.h"



// cost1 ---------
// NL16 ------------
#include "model_include/IK_residual_func_cost1_Nl16.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl16.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl16.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl16.h"


// cost2 ---------
// Nl16 ------------
#include "model_include/IK_residual_func_cost2_Nl16.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl16.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl16.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl16.h"


// cost3 ---------
// Nl16 ------------
#include "model_include/IK_residual_func_cost3_Nl16.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl16.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl16.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl16.h"



// cost1 ---------
// NL17 ------------
#include "model_include/IK_residual_func_cost1_Nl17.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl17.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl17.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl17.h"


// cost2 ---------
// Nl17 ------------
#include "model_include/IK_residual_func_cost2_Nl17.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl17.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl17.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl17.h"


// cost3 ---------
// Nl17 ------------
#include "model_include/IK_residual_func_cost3_Nl17.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl17.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl17.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl17.h"



// cost1 ---------
// NL18 ------------
#include "model_include/IK_residual_func_cost1_Nl18.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl18.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl18.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl18.h"


// cost2 ---------
// Nl18 ------------
#include "model_include/IK_residual_func_cost2_Nl18.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl18.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl18.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl18.h"


// cost3 ---------
// Nl18 ------------
#include "model_include/IK_residual_func_cost3_Nl18.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl18.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl18.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl18.h"



// cost1 ---------
// NL19 ------------
#include "model_include/IK_residual_func_cost1_Nl19.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl19.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl19.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl19.h"


// cost2 ---------
// Nl19 ------------
#include "model_include/IK_residual_func_cost2_Nl19.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl19.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl19.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl19.h"


// cost3 ---------
// Nl19 ------------
#include "model_include/IK_residual_func_cost3_Nl19.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl19.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl19.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl19.h"



// cost1 ---------
// NL20 ------------
#include "model_include/IK_residual_func_cost1_Nl20.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl20.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl20.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl20.h"


// cost2 ---------
// Nl20 ------------
#include "model_include/IK_residual_func_cost2_Nl20.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl20.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl20.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl20.h"


// cost3 ---------
// Nl20 ------------
#include "model_include/IK_residual_func_cost3_Nl20.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl20.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl20.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl20.h"



// cost1 ---------
// NL21 ------------
#include "model_include/IK_residual_func_cost1_Nl21.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl21.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl21.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl21.h"


// cost2 ---------
// Nl21 ------------
#include "model_include/IK_residual_func_cost2_Nl21.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl21.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl21.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl21.h"


// cost3 ---------
// Nl21 ------------
#include "model_include/IK_residual_func_cost3_Nl21.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl21.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl21.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl21.h"



// cost1 ---------
// NL22 ------------
#include "model_include/IK_residual_func_cost1_Nl22.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl22.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl22.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl22.h"


// cost2 ---------
// Nl22 ------------
#include "model_include/IK_residual_func_cost2_Nl22.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl22.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl22.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl22.h"


// cost3 ---------
// Nl22 ------------
#include "model_include/IK_residual_func_cost3_Nl22.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl22.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl22.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl22.h"



// cost1 ---------
// NL23 ------------
#include "model_include/IK_residual_func_cost1_Nl23.h"
#include "model_include/IK_residual_func_cost1_wrt_fh1_Nl23.h"
#include "model_include/IK_residual_func_cost1_wrt_fv1_Nl23.h"
#include "model_include/IK_residual_func_cost1_wrt_DeltaRot_Nl23.h"


// cost2 ---------
// Nl23 ------------
#include "model_include/IK_residual_func_cost2_Nl23.h"
#include "model_include/IK_residual_func_cost2_wrt_fh1_Nl23.h"
#include "model_include/IK_residual_func_cost2_wrt_fv1_Nl23.h"
#include "model_include/IK_residual_func_cost2_wrt_DeltaRot_Nl23.h"


// cost3 ---------
// Nl23 ------------
#include "model_include/IK_residual_func_cost3_Nl23.h"
#include "model_include/IK_residual_func_cost3_wrt_fh1_Nl23.h"
#include "model_include/IK_residual_func_cost3_wrt_fv1_Nl23.h"
#include "model_include/IK_residual_func_cost3_wrt_DeltaRot_Nl23.h"


