// -----------------------------------------------------------------------------
// This file was autogenerated by symforce from template:
//     function/FUNCTION.h.jinja
// Do NOT modify by hand.
// -----------------------------------------------------------------------------

#pragma once

#include <Eigen/Dense>

#include <sym/rot3.h>

namespace sym {

/**
 * This function was autogenerated from a symbolic function. Do not modify by hand.
 *
 * Symbolic function: IK_residual_func_cost3_Nl12
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     position_vector: Matrix31
 *     encoder: Matrix41
 *     p_a: Matrix31
 *     p_b: Matrix31
 *     p_c: Matrix31
 *     p_d: Matrix31
 *     Rot_init: Rot3
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix41
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost3Nl12(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 285

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (104)
  const Scalar _tmp0 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp2 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp3 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp4 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp4 +
                       Scalar(0.20999999999999999);
  const Scalar _tmp6 = -_tmp5;
  const Scalar _tmp7 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp8 = 2 * _tmp3 * _tmp7;
  const Scalar _tmp9 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp10 = _tmp1 * _tmp9;
  const Scalar _tmp11 =
      -Scalar(0.010999999999999999) * _tmp10 - Scalar(0.010999999999999999) * _tmp8;
  const Scalar _tmp12 = 2 * _tmp1;
  const Scalar _tmp13 = _tmp12 * _tmp7;
  const Scalar _tmp14 = _tmp3 * _tmp9;
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 - Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp11 + _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp6;
  const Scalar _tmp18 = _tmp17 - p_d(0, 0) + position_vector(0, 0);
  const Scalar _tmp19 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp20 = -_tmp19;
  const Scalar _tmp21 = 1 - 2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp21 + Scalar(0.20999999999999999) * _tmp4;
  const Scalar _tmp23 = _tmp12 * _tmp3;
  const Scalar _tmp24 = _tmp7 * _tmp9;
  const Scalar _tmp25 =
      -Scalar(0.010999999999999999) * _tmp23 + Scalar(0.010999999999999999) * _tmp24;
  const Scalar _tmp26 = _tmp22 + _tmp25;
  const Scalar _tmp27 = _tmp20 + _tmp26;
  const Scalar _tmp28 = _tmp27 - p_d(1, 0) + position_vector(1, 0);
  const Scalar _tmp29 = std::pow(Scalar(std::pow(_tmp18, Scalar(2)) + std::pow(_tmp28, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp30 = _tmp18 * _tmp29;
  const Scalar _tmp31 = -Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp32 = -_tmp31;
  const Scalar _tmp33 =
      -Scalar(0.010999999999999999) * _tmp2 - Scalar(0.010999999999999999) * _tmp21;
  const Scalar _tmp34 = Scalar(0.20999999999999999) * _tmp23 + Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp35 = _tmp33 - _tmp34;
  const Scalar _tmp36 = _tmp32 + _tmp35;
  const Scalar _tmp37 = _tmp30 * _tmp36;
  const Scalar _tmp38 = _tmp31 + _tmp35;
  const Scalar _tmp39 = _tmp11 - _tmp15;
  const Scalar _tmp40 = _tmp39 + _tmp5;
  const Scalar _tmp41 = _tmp40 - p_b(0, 0) + position_vector(0, 0);
  const Scalar _tmp42 = -_tmp22 + _tmp25;
  const Scalar _tmp43 = _tmp19 + _tmp42;
  const Scalar _tmp44 = _tmp43 - p_b(1, 0) + position_vector(1, 0);
  const Scalar _tmp45 = std::pow(Scalar(std::pow(_tmp41, Scalar(2)) + std::pow(_tmp44, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp46 = _tmp41 * _tmp45;
  const Scalar _tmp47 = _tmp36 * _tmp46;
  const Scalar _tmp48 = -_tmp38 * _tmp46 + _tmp47;
  const Scalar _tmp49 = _tmp44 * _tmp45;
  const Scalar _tmp50 = _tmp39 + _tmp6;
  const Scalar _tmp51 = _tmp50 - p_a(0, 0) + position_vector(0, 0);
  const Scalar _tmp52 = Scalar(1.0) / (_tmp51);
  const Scalar _tmp53 = _tmp20 + _tmp42;
  const Scalar _tmp54 = _tmp53 - p_a(1, 0) + position_vector(1, 0);
  const Scalar _tmp55 = _tmp52 * _tmp54;
  const Scalar _tmp56 = Scalar(1.0) / (_tmp46 * _tmp55 - _tmp49);
  const Scalar _tmp57 = _tmp28 * _tmp29;
  const Scalar _tmp58 = _tmp30 * _tmp55 - _tmp57;
  const Scalar _tmp59 = _tmp56 * _tmp58;
  const Scalar _tmp60 = _tmp32 + _tmp33 + _tmp34;
  const Scalar _tmp61 = _tmp38 * _tmp49 - _tmp47 * _tmp55;
  const Scalar _tmp62 = Scalar(1.0) * _tmp53;
  const Scalar _tmp63 = Scalar(1.0) * _tmp50;
  const Scalar _tmp64 = (-_tmp40 + _tmp63) / (_tmp43 - _tmp62);
  const Scalar _tmp65 = -_tmp30 * _tmp60 + _tmp37 - _tmp48 * _tmp59 -
                        _tmp64 * (-_tmp37 * _tmp55 + _tmp57 * _tmp60 - _tmp59 * _tmp61);
  const Scalar _tmp66 = Scalar(1.0) / (_tmp65);
  const Scalar _tmp67 = _tmp62 * _tmp64 + _tmp63;
  const Scalar _tmp68 = 0;
  const Scalar _tmp69 = _tmp46 * _tmp59;
  const Scalar _tmp70 =
      std::sqrt(Scalar(std::pow(_tmp51, Scalar(2)) + std::pow(_tmp54, Scalar(2))));
  const Scalar _tmp71 = _tmp52 * _tmp70;
  const Scalar _tmp72 = Scalar(1.0) / (_tmp70);
  const Scalar _tmp73 = _tmp71 * (_tmp50 * _tmp54 * _tmp72 - _tmp51 * _tmp53 * _tmp72);
  const Scalar _tmp74 = -_tmp40 * _tmp49 + _tmp43 * _tmp46 + _tmp46 * _tmp73;
  const Scalar _tmp75 = _tmp55 * _tmp56;
  const Scalar _tmp76 = -_tmp36 + _tmp48 * _tmp75 - _tmp64 * (_tmp36 * _tmp55 + _tmp61 * _tmp75);
  const Scalar _tmp77 = -_tmp17 * _tmp57 + _tmp27 * _tmp30 + _tmp30 * _tmp73 - _tmp59 * _tmp74;
  const Scalar _tmp78 = _tmp66 * _tmp77;
  const Scalar _tmp79 = Scalar(1.0) / (_tmp77);
  const Scalar _tmp80 = _tmp65 * _tmp79;
  const Scalar _tmp81 = _tmp76 + _tmp80 * (-_tmp73 + _tmp74 * _tmp75 - _tmp76 * _tmp78);
  const Scalar _tmp82 = _tmp58 * _tmp66;
  const Scalar _tmp83 = -_tmp55 - _tmp81 * _tmp82;
  const Scalar _tmp84 = _tmp46 * _tmp56;
  const Scalar _tmp85 = _tmp30 * _tmp66;
  const Scalar _tmp86 = _tmp16 + _tmp5;
  const Scalar _tmp87 = _tmp86 - p_c(0, 0) + position_vector(0, 0);
  const Scalar _tmp88 = _tmp19 + _tmp26;
  const Scalar _tmp89 = _tmp88 - p_c(1, 0) + position_vector(1, 0);
  const Scalar _tmp90 = std::pow(Scalar(std::pow(_tmp87, Scalar(2)) + std::pow(_tmp89, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp91 = _tmp87 * _tmp90;
  const Scalar _tmp92 = _tmp71 * fh1;
  const Scalar _tmp93 = Scalar(1.0) * _tmp79;
  const Scalar _tmp94 = _tmp89 * _tmp90;
  const Scalar _tmp95 = fh1 * (_tmp86 * _tmp94 - _tmp88 * _tmp91);
  const Scalar _tmp96 = Scalar(1.0) * _tmp56;
  const Scalar _tmp97 = -_tmp48 * _tmp96 + _tmp61 * _tmp64 * _tmp96;
  const Scalar _tmp98 = _tmp80 * (-_tmp74 * _tmp96 - _tmp78 * _tmp97) + _tmp97;
  const Scalar _tmp99 = -_tmp82 * _tmp98 + Scalar(1.0);
  const Scalar _tmp100 = _tmp56 * fh1;
  const Scalar _tmp101 = _tmp0 * _tmp68;
  const Scalar _tmp102 = _tmp93 * _tmp95;
  const Scalar _tmp103 = _tmp66 * fh1;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = std::exp(-fh1);
  _res(1, 0) = std::exp(_tmp0 * _tmp71 * (_tmp30 * _tmp68 - _tmp68 * _tmp69) +
                        _tmp71 * _tmp95 * (_tmp30 * _tmp93 - _tmp69 * _tmp93) +
                        _tmp91 * _tmp92 * (_tmp81 * _tmp85 + _tmp83 * _tmp84 + Scalar(1.0)) +
                        _tmp92 * _tmp94 * (_tmp84 * _tmp99 + _tmp85 * _tmp98));
  _res(2, 0) = std::exp(-_tmp100 * _tmp83 * _tmp91 - _tmp100 * _tmp94 * _tmp99 + _tmp101 * _tmp59 +
                        _tmp102 * _tmp59);
  _res(3, 0) = std::exp(-_tmp101 - _tmp102 - _tmp103 * _tmp81 * _tmp91 - _tmp103 * _tmp94 * _tmp98);

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
