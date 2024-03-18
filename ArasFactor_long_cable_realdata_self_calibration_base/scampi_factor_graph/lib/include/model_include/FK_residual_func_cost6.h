// -----------------------------------------------------------------------------
// This file was autogenerated by symforce from template:
//     function/FUNCTION.h.jinja
// Do NOT modify by hand.
// -----------------------------------------------------------------------------

#pragma once

#include <Eigen/Dense>

#include <sym/pose3.h>
#include <sym/rot3.h>

namespace sym {

/**
 * This function was autogenerated from a symbolic function. Do not modify by hand.
 *
 * Symbolic function: FK_residual_func_cost6
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     TransformationMatrix: Pose3
 *     encoder: Matrix41
 *     p_a: Matrix31
 *     p_b: Matrix31
 *     p_c: Matrix31
 *     p_d: Matrix31
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix21
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 2, 1> FkResidualFuncCost6(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const sym::Pose3<Scalar>& TransformationMatrix, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const Scalar epsilon) {
  // Total ops: 304

  // Unused inputs
  (void)fh1;
  (void)fv1;
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 7, 1>& _TransformationMatrix = TransformationMatrix.Data();

  // Intermediate terms (110)
  const Scalar _tmp0 =
      _DeltaRot[0] * _TransformationMatrix[3] - _DeltaRot[1] * _TransformationMatrix[2] +
      _DeltaRot[2] * _TransformationMatrix[1] + _DeltaRot[3] * _TransformationMatrix[0];
  const Scalar _tmp1 =
      -_DeltaRot[0] * _TransformationMatrix[1] + _DeltaRot[1] * _TransformationMatrix[0] +
      _DeltaRot[2] * _TransformationMatrix[3] + _DeltaRot[3] * _TransformationMatrix[2];
  const Scalar _tmp2 = 2 * _tmp1;
  const Scalar _tmp3 = _tmp0 * _tmp2;
  const Scalar _tmp4 =
      _DeltaRot[0] * _TransformationMatrix[2] + _DeltaRot[1] * _TransformationMatrix[3] -
      _DeltaRot[2] * _TransformationMatrix[0] + _DeltaRot[3] * _TransformationMatrix[1];
  const Scalar _tmp5 =
      -2 * _DeltaRot[0] * _TransformationMatrix[0] - 2 * _DeltaRot[1] * _TransformationMatrix[1] -
      2 * _DeltaRot[2] * _TransformationMatrix[2] + 2 * _DeltaRot[3] * _TransformationMatrix[3];
  const Scalar _tmp6 = _tmp4 * _tmp5;
  const Scalar _tmp7 = _tmp3 + _tmp6;
  const Scalar _tmp8 = Scalar(0.20999999999999999) * _tmp3 - Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp9 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp10 = 1 - 2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp11 =
      -Scalar(0.010999999999999999) * _tmp10 - Scalar(0.010999999999999999) * _tmp9;
  const Scalar _tmp12 = _tmp2 * _tmp4;
  const Scalar _tmp13 = _tmp0 * _tmp5;
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp11 - _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp8;
  const Scalar _tmp17 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp18 = Scalar(0.20999999999999999) * _tmp17 + Scalar(0.20999999999999999) * _tmp9 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp19 = -_tmp18;
  const Scalar _tmp20 = 2 * _tmp0 * _tmp4;
  const Scalar _tmp21 = _tmp1 * _tmp5;
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp20 + Scalar(0.20999999999999999) * _tmp21;
  const Scalar _tmp23 = _tmp12 - _tmp13;
  const Scalar _tmp24 = -Scalar(0.010999999999999999) * _tmp23;
  const Scalar _tmp25 = _tmp22 + _tmp24;
  const Scalar _tmp26 = _tmp19 + _tmp25;
  const Scalar _tmp27 = _TransformationMatrix[5] + _tmp26 - p_b(1, 0);
  const Scalar _tmp28 = _TransformationMatrix[6] + _tmp16 - p_b(2, 0);
  const Scalar _tmp29 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp17;
  const Scalar _tmp30 = -Scalar(0.010999999999999999) * _tmp7;
  const Scalar _tmp31 = Scalar(0.20999999999999999) * _tmp20 - Scalar(0.20999999999999999) * _tmp21;
  const Scalar _tmp32 = _tmp30 - _tmp31;
  const Scalar _tmp33 = _tmp29 + _tmp32;
  const Scalar _tmp34 = _TransformationMatrix[4] + _tmp33 - p_b(0, 0);
  const Scalar _tmp35 = std::pow(Scalar(std::pow(_tmp27, Scalar(2)) + std::pow(_tmp28, Scalar(2)) +
                                        std::pow(_tmp34, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp36 = _tmp34 * _tmp35;
  const Scalar _tmp37 = _tmp28 * _tmp35;
  const Scalar _tmp38 = _tmp16 * _tmp36 - _tmp33 * _tmp37;
  const Scalar _tmp39 = -_tmp22 + _tmp24;
  const Scalar _tmp40 = _tmp19 + _tmp39;
  const Scalar _tmp41 = _TransformationMatrix[5] + _tmp40 - p_a(1, 0);
  const Scalar _tmp42 = -_tmp29;
  const Scalar _tmp43 = _tmp32 + _tmp42;
  const Scalar _tmp44 = _TransformationMatrix[4] + _tmp43 - p_a(0, 0);
  const Scalar _tmp45 = Scalar(1.0) / (_tmp44);
  const Scalar _tmp46 = _tmp36 * _tmp45;
  const Scalar _tmp47 = _tmp27 * _tmp35;
  const Scalar _tmp48 = Scalar(1.0) / (-_tmp41 * _tmp46 + _tmp47);
  const Scalar _tmp49 = _tmp11 + _tmp14;
  const Scalar _tmp50 = _tmp49 + _tmp8;
  const Scalar _tmp51 = _TransformationMatrix[6] + _tmp50 - p_c(2, 0);
  const Scalar _tmp52 = _tmp30 + _tmp31;
  const Scalar _tmp53 = _tmp29 + _tmp52;
  const Scalar _tmp54 = _TransformationMatrix[4] + _tmp53 - p_c(0, 0);
  const Scalar _tmp55 = _tmp18 + _tmp25;
  const Scalar _tmp56 = _TransformationMatrix[5] + _tmp55 - p_c(1, 0);
  const Scalar _tmp57 = std::pow(Scalar(std::pow(_tmp51, Scalar(2)) + std::pow(_tmp54, Scalar(2)) +
                                        std::pow(_tmp56, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp58 = _tmp54 * _tmp57;
  const Scalar _tmp59 = _tmp41 * _tmp45;
  const Scalar _tmp60 = _tmp56 * _tmp57;
  const Scalar _tmp61 = -_tmp58 * _tmp59 + _tmp60;
  const Scalar _tmp62 = -_tmp8;
  const Scalar _tmp63 = _tmp49 + _tmp62;
  const Scalar _tmp64 = _tmp18 + _tmp39;
  const Scalar _tmp65 = _TransformationMatrix[5] + _tmp64 - p_d(1, 0);
  const Scalar _tmp66 = _TransformationMatrix[6] + _tmp63 - p_d(2, 0);
  const Scalar _tmp67 = _tmp42 + _tmp52;
  const Scalar _tmp68 = _TransformationMatrix[4] + _tmp67 - p_d(0, 0);
  const Scalar _tmp69 = std::pow(Scalar(std::pow(_tmp65, Scalar(2)) + std::pow(_tmp66, Scalar(2)) +
                                        std::pow(_tmp68, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp70 = _tmp65 * _tmp69;
  const Scalar _tmp71 = _tmp68 * _tmp69;
  const Scalar _tmp72 = _tmp15 + _tmp62;
  const Scalar _tmp73 = _TransformationMatrix[6] + _tmp72 - p_a(2, 0);
  const Scalar _tmp74 = _tmp45 * _tmp73;
  const Scalar _tmp75 = _tmp66 * _tmp69;
  const Scalar _tmp76 = _tmp37 - _tmp46 * _tmp73;
  const Scalar _tmp77 = -_tmp59 * _tmp71 + _tmp70;
  const Scalar _tmp78 = _tmp48 * _tmp77;
  const Scalar _tmp79 = -_tmp71 * _tmp74 + _tmp75 - _tmp76 * _tmp78;
  const Scalar _tmp80 = _tmp51 * _tmp57;
  const Scalar _tmp81 = _tmp48 * _tmp61;
  const Scalar _tmp82 = Scalar(1.0) / (-_tmp58 * _tmp74 - _tmp76 * _tmp81 + _tmp80);
  const Scalar _tmp83 = std::sqrt(Scalar(std::pow(_tmp41, Scalar(2)) + std::pow(_tmp44, Scalar(2)) +
                                         std::pow(_tmp73, Scalar(2))));
  const Scalar _tmp84 = Scalar(1.0) / (_tmp83);
  const Scalar _tmp85 = _tmp73 * _tmp84;
  const Scalar _tmp86 = _tmp41 * _tmp84;
  const Scalar _tmp87 = _tmp45 * _tmp83;
  const Scalar _tmp88 = _tmp87 * (_tmp40 * _tmp85 - _tmp72 * _tmp86);
  const Scalar _tmp89 = -_tmp16 * _tmp47 + _tmp26 * _tmp37 - _tmp36 * _tmp88;
  const Scalar _tmp90 =
      _tmp82 * (-_tmp50 * _tmp60 + _tmp55 * _tmp80 - _tmp58 * _tmp88 - _tmp81 * _tmp89);
  const Scalar _tmp91 = _tmp79 * _tmp90;
  const Scalar _tmp92 = Scalar(1.0) / (-_tmp63 * _tmp70 + _tmp64 * _tmp75 - _tmp71 * _tmp88 -
                                       _tmp78 * _tmp89 - _tmp91);
  const Scalar _tmp93 = _tmp79 * _tmp82 * _tmp92;
  const Scalar _tmp94 = _tmp77 * _tmp92;
  const Scalar _tmp95 = _tmp48 * (_tmp61 * _tmp93 - _tmp94);
  const Scalar _tmp96 = _tmp92 * (_tmp63 * _tmp71 - _tmp67 * _tmp75);
  const Scalar _tmp97 = _tmp71 * _tmp92;
  const Scalar _tmp98 = -_tmp36 * _tmp95 + _tmp58 * _tmp93 - _tmp97;
  const Scalar _tmp99 = _tmp44 * _tmp84;
  const Scalar _tmp100 = _tmp87 * (-_tmp43 * _tmp85 + _tmp72 * _tmp99);
  const Scalar _tmp101 = _tmp50 * _tmp58 - _tmp53 * _tmp80;
  const Scalar _tmp102 = Scalar(40.024799999999999) * _tmp23;
  const Scalar _tmp103 = _tmp82 * (_tmp91 * _tmp92 + 1);
  const Scalar _tmp104 = _tmp48 * (-_tmp103 * _tmp61 + _tmp90 * _tmp94);
  const Scalar _tmp105 = -_tmp103 * _tmp58 - _tmp104 * _tmp36 + _tmp90 * _tmp97;
  const Scalar _tmp106 = _tmp53 * _tmp60 - _tmp55 * _tmp58;
  const Scalar _tmp107 = _tmp92 * (-_tmp64 * _tmp71 + _tmp67 * _tmp70);
  const Scalar _tmp108 = -_tmp26 * _tmp36 + _tmp33 * _tmp47;
  const Scalar _tmp109 = _tmp87 * (-_tmp40 * _tmp99 + _tmp43 * _tmp86);

  // Output terms (1)
  Eigen::Matrix<Scalar, 2, 1> _res;

  _res(0, 0) = Scalar(333.54000000000002) * _tmp100 * _tmp105 +
               Scalar(333.54000000000002) * _tmp101 * _tmp103 -
               _tmp102 * (_tmp100 * _tmp98 - _tmp101 * _tmp93 + _tmp38 * _tmp95 + _tmp96) +
               Scalar(333.54000000000002) * _tmp104 * _tmp38 - Scalar(40.024799999999999) * _tmp7 -
               Scalar(333.54000000000002) * _tmp90 * _tmp96;
  _res(1, 0) = -_tmp102 * (-_tmp106 * _tmp93 + _tmp107 + _tmp108 * _tmp95 + _tmp109 * _tmp98) +
               Scalar(333.54000000000002) * _tmp103 * _tmp106 +
               Scalar(333.54000000000002) * _tmp104 * _tmp108 +
               Scalar(333.54000000000002) * _tmp105 * _tmp109 -
               Scalar(333.54000000000002) * _tmp107 * _tmp90;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
