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
 * Symbolic function: IK_residual_func_cost3_wrt_fh1_Nl18
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost3WrtFh1Nl18(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 305

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (108)
  const Scalar _tmp0 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp2 = 2 * _tmp1;
  const Scalar _tmp3 = _tmp0 * _tmp2;
  const Scalar _tmp4 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp5 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp6 = _tmp4 * _tmp5;
  const Scalar _tmp7 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp8 = -_tmp7;
  const Scalar _tmp9 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp10 = 1 - 2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp11 =
      -Scalar(0.010999999999999999) * _tmp10 - Scalar(0.010999999999999999) * _tmp9;
  const Scalar _tmp12 = 2 * _tmp0 * _tmp4;
  const Scalar _tmp13 = _tmp1 * _tmp5;
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 - Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp11 + _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp8;
  const Scalar _tmp17 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp18 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp17;
  const Scalar _tmp19 =
      -Scalar(0.010999999999999999) * _tmp12 - Scalar(0.010999999999999999) * _tmp13;
  const Scalar _tmp20 = _tmp2 * _tmp4;
  const Scalar _tmp21 = _tmp0 * _tmp5;
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp20 - Scalar(0.20999999999999999) * _tmp21;
  const Scalar _tmp23 = _tmp19 - _tmp22;
  const Scalar _tmp24 = _tmp18 + _tmp23;
  const Scalar _tmp25 = _tmp24 - p_b(0, 0) + position_vector(0, 0);
  const Scalar _tmp26 = Scalar(0.20999999999999999) * _tmp20 + Scalar(0.20999999999999999) * _tmp21;
  const Scalar _tmp27 = Scalar(0.20999999999999999) * _tmp17 + Scalar(0.20999999999999999) * _tmp9 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp28 =
      -Scalar(0.010999999999999999) * _tmp3 + Scalar(0.010999999999999999) * _tmp6;
  const Scalar _tmp29 = -_tmp27 + _tmp28;
  const Scalar _tmp30 = _tmp26 + _tmp29;
  const Scalar _tmp31 = _tmp30 - p_b(1, 0) + position_vector(1, 0);
  const Scalar _tmp32 = std::pow(Scalar(std::pow(_tmp25, Scalar(2)) + std::pow(_tmp31, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp33 = _tmp31 * _tmp32;
  const Scalar _tmp34 = -_tmp18;
  const Scalar _tmp35 = _tmp23 + _tmp34;
  const Scalar _tmp36 = _tmp35 - p_a(0, 0) + position_vector(0, 0);
  const Scalar _tmp37 = Scalar(1.0) / (_tmp36);
  const Scalar _tmp38 = -_tmp26;
  const Scalar _tmp39 = _tmp29 + _tmp38;
  const Scalar _tmp40 = _tmp39 - p_a(1, 0) + position_vector(1, 0);
  const Scalar _tmp41 = _tmp37 * _tmp40;
  const Scalar _tmp42 = _tmp11 - _tmp14 + _tmp8;
  const Scalar _tmp43 = _tmp25 * _tmp32;
  const Scalar _tmp44 = _tmp42 * _tmp43;
  const Scalar _tmp45 = _tmp16 * _tmp33 - _tmp41 * _tmp44;
  const Scalar _tmp46 = Scalar(1.0) * _tmp39;
  const Scalar _tmp47 = Scalar(1.0) * _tmp35;
  const Scalar _tmp48 = (-_tmp24 + _tmp47) / (_tmp30 - _tmp46);
  const Scalar _tmp49 = Scalar(1.0) / (-_tmp33 + _tmp41 * _tmp43);
  const Scalar _tmp50 = Scalar(1.0) * _tmp49;
  const Scalar _tmp51 = -_tmp16 * _tmp43 + _tmp44;
  const Scalar _tmp52 = _tmp45 * _tmp48 * _tmp50 - _tmp50 * _tmp51;
  const Scalar _tmp53 = _tmp27 + _tmp28;
  const Scalar _tmp54 = _tmp26 + _tmp53;
  const Scalar _tmp55 = _tmp54 - p_c(1, 0) + position_vector(1, 0);
  const Scalar _tmp56 = _tmp19 + _tmp22;
  const Scalar _tmp57 = _tmp18 + _tmp56;
  const Scalar _tmp58 = _tmp57 - p_c(0, 0) + position_vector(0, 0);
  const Scalar _tmp59 = std::pow(Scalar(std::pow(_tmp55, Scalar(2)) + std::pow(_tmp58, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp60 = _tmp55 * _tmp59;
  const Scalar _tmp61 = _tmp58 * _tmp59;
  const Scalar _tmp62 = _tmp41 * _tmp61 - _tmp60;
  const Scalar _tmp63 = _tmp49 * _tmp62;
  const Scalar _tmp64 = _tmp42 * _tmp61;
  const Scalar _tmp65 = _tmp15 + _tmp7;
  const Scalar _tmp66 = -_tmp48 * (-_tmp41 * _tmp64 - _tmp45 * _tmp63 + _tmp60 * _tmp65) -
                        _tmp51 * _tmp63 - _tmp61 * _tmp65 + _tmp64;
  const Scalar _tmp67 = Scalar(1.0) / (_tmp66);
  const Scalar _tmp68 =
      std::sqrt(Scalar(std::pow(_tmp36, Scalar(2)) + std::pow(_tmp40, Scalar(2))));
  const Scalar _tmp69 = Scalar(1.0) / (_tmp68);
  const Scalar _tmp70 = _tmp37 * _tmp68;
  const Scalar _tmp71 = _tmp70 * (_tmp35 * _tmp40 * _tmp69 - _tmp36 * _tmp39 * _tmp69);
  const Scalar _tmp72 = -_tmp24 * _tmp33 + _tmp30 * _tmp43 + _tmp43 * _tmp71;
  const Scalar _tmp73 = _tmp54 * _tmp61 - _tmp57 * _tmp60 + _tmp61 * _tmp71 - _tmp63 * _tmp72;
  const Scalar _tmp74 = _tmp67 * _tmp73;
  const Scalar _tmp75 = Scalar(1.0) / (_tmp73);
  const Scalar _tmp76 = _tmp66 * _tmp75;
  const Scalar _tmp77 = _tmp52 + _tmp76 * (-_tmp50 * _tmp72 - _tmp52 * _tmp74);
  const Scalar _tmp78 = _tmp61 * _tmp67;
  const Scalar _tmp79 = _tmp62 * _tmp67;
  const Scalar _tmp80 = _tmp49 * (-_tmp77 * _tmp79 + Scalar(1.0));
  const Scalar _tmp81 = _tmp34 + _tmp56;
  const Scalar _tmp82 = _tmp81 - p_d(0, 0) + position_vector(0, 0);
  const Scalar _tmp83 = _tmp38 + _tmp53;
  const Scalar _tmp84 = _tmp83 - p_d(1, 0) + position_vector(1, 0);
  const Scalar _tmp85 = std::pow(Scalar(std::pow(_tmp82, Scalar(2)) + std::pow(_tmp84, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp86 = _tmp84 * _tmp85;
  const Scalar _tmp87 = _tmp70 * _tmp86 * (_tmp43 * _tmp80 + _tmp77 * _tmp78);
  const Scalar _tmp88 = _tmp82 * _tmp85;
  const Scalar _tmp89 = _tmp81 * _tmp86 - _tmp83 * _tmp88;
  const Scalar _tmp90 = _tmp50 * _tmp62 * _tmp75;
  const Scalar _tmp91 = Scalar(1.0) * _tmp75;
  const Scalar _tmp92 = _tmp70 * _tmp89 * (-_tmp43 * _tmp90 + _tmp61 * _tmp91);
  const Scalar _tmp93 = _tmp41 * _tmp49;
  const Scalar _tmp94 = -_tmp42 - _tmp48 * (_tmp41 * _tmp42 + _tmp45 * _tmp93) + _tmp51 * _tmp93;
  const Scalar _tmp95 = _tmp76 * (-_tmp71 + _tmp72 * _tmp93 - _tmp74 * _tmp94) + _tmp94;
  const Scalar _tmp96 = _tmp49 * (-_tmp41 - _tmp79 * _tmp95);
  const Scalar _tmp97 = _tmp70 * _tmp88 * (_tmp43 * _tmp96 + _tmp78 * _tmp95 + Scalar(1.0));
  const Scalar _tmp98 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp99 = _tmp46 * _tmp48 + _tmp47;
  const Scalar _tmp100 = 0;
  const Scalar _tmp101 = _tmp88 * _tmp96;
  const Scalar _tmp102 = _tmp80 * _tmp86;
  const Scalar _tmp103 = _tmp89 * _tmp90;
  const Scalar _tmp104 = _tmp100 * _tmp98;
  const Scalar _tmp105 = _tmp89 * _tmp91;
  const Scalar _tmp106 = _tmp67 * _tmp88 * _tmp95;
  const Scalar _tmp107 = _tmp67 * _tmp77 * _tmp86;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = -std::exp(-fh1);
  _res(1, 0) = -(-_tmp87 - _tmp92 - _tmp97) *
               std::exp(_tmp70 * _tmp98 * (-_tmp100 * _tmp43 * _tmp63 + _tmp100 * _tmp61) +
                        _tmp87 * fh1 + _tmp92 * fh1 + _tmp97 * fh1);
  _res(2, 0) = -(_tmp101 + _tmp102 - _tmp103) *
               std::exp(-_tmp101 * fh1 - _tmp102 * fh1 + _tmp103 * fh1 + _tmp104 * _tmp63);
  _res(3, 0) = -(_tmp105 + _tmp106 + _tmp107) *
               std::exp(-_tmp104 - _tmp105 * fh1 - _tmp106 * fh1 - _tmp107 * fh1);

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
