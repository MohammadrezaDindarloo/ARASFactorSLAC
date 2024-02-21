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
 * Symbolic function: IK_residual_func_cost3_wrt_fv1_Nl11
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     position_vector: Matrix31
 *     Rot_init: Rot3
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix41
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost3WrtFv1Nl11(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 281

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (109)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp1 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp2 = 2 * _tmp1;
  const Scalar _tmp3 = _tmp0 * _tmp2;
  const Scalar _tmp4 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp5 = -_DeltaRot[0] * _Rot_init[0] - _DeltaRot[1] * _Rot_init[1] -
                       _DeltaRot[2] * _Rot_init[2] + _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp6 = 2 * _tmp4 * _tmp5;
  const Scalar _tmp7 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp8 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp9 = 1 - 2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp10 =
      -Scalar(0.010999999999999999) * _tmp8 - Scalar(0.010999999999999999) * _tmp9;
  const Scalar _tmp11 = _tmp2 * _tmp4;
  const Scalar _tmp12 = 2 * _tmp0;
  const Scalar _tmp13 = _tmp12 * _tmp5;
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp11 - Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp10 - _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp7;
  const Scalar _tmp17 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp18 = Scalar(0.20999999999999999) * _tmp17 + Scalar(0.20999999999999999) * _tmp8 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp19 =
      -Scalar(0.010999999999999999) * _tmp3 + Scalar(0.010999999999999999) * _tmp6;
  const Scalar _tmp20 = _tmp12 * _tmp4;
  const Scalar _tmp21 = _tmp2 * _tmp5;
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp20 + Scalar(0.20999999999999999) * _tmp21;
  const Scalar _tmp23 = _tmp19 - _tmp22;
  const Scalar _tmp24 = _tmp18 + _tmp23;
  const Scalar _tmp25 = position_vector(1, 0) + Scalar(-110.0);
  const Scalar _tmp26 = _tmp24 + _tmp25;
  const Scalar _tmp27 = Scalar(0.20999999999999999) * _tmp20 - Scalar(0.20999999999999999) * _tmp21;
  const Scalar _tmp28 = Scalar(0.20999999999999999) * _tmp17 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp29 =
      -Scalar(0.010999999999999999) * _tmp11 - Scalar(0.010999999999999999) * _tmp13;
  const Scalar _tmp30 = -_tmp28 + _tmp29;
  const Scalar _tmp31 = _tmp27 + _tmp30;
  const Scalar _tmp32 = position_vector(0, 0) + Scalar(125.0);
  const Scalar _tmp33 = _tmp31 + _tmp32;
  const Scalar _tmp34 = Scalar(1.0) / (_tmp33);
  const Scalar _tmp35 = _tmp26 * _tmp34;
  const Scalar _tmp36 = _tmp16 * _tmp35;
  const Scalar _tmp37 = _tmp19 + _tmp22;
  const Scalar _tmp38 = _tmp18 + _tmp37;
  const Scalar _tmp39 = _tmp25 + _tmp38;
  const Scalar _tmp40 = _tmp28 + _tmp29;
  const Scalar _tmp41 = _tmp27 + _tmp40;
  const Scalar _tmp42 = position_vector(0, 0) + Scalar(-125.0);
  const Scalar _tmp43 = _tmp41 + _tmp42;
  const Scalar _tmp44 = std::pow(Scalar(std::pow(_tmp39, Scalar(2)) + std::pow(_tmp43, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp45 = _tmp43 * _tmp44;
  const Scalar _tmp46 = _tmp16 * _tmp45;
  const Scalar _tmp47 = _tmp10 + _tmp14 + _tmp7;
  const Scalar _tmp48 = _tmp39 * _tmp44;
  const Scalar _tmp49 = -_tmp35 * _tmp46 + _tmp47 * _tmp48;
  const Scalar _tmp50 = Scalar(1.0) / (_tmp35 * _tmp45 - _tmp48);
  const Scalar _tmp51 = _tmp35 * _tmp50;
  const Scalar _tmp52 = Scalar(1.0) * _tmp24;
  const Scalar _tmp53 = Scalar(1.0) * _tmp31;
  const Scalar _tmp54 = (-_tmp41 + _tmp53) / (_tmp38 - _tmp52);
  const Scalar _tmp55 = -_tmp45 * _tmp47 + _tmp46;
  const Scalar _tmp56 = -_tmp16 + _tmp51 * _tmp55 - _tmp54 * (_tmp36 + _tmp49 * _tmp51);
  const Scalar _tmp57 = -_tmp18;
  const Scalar _tmp58 = _tmp23 + _tmp57;
  const Scalar _tmp59 = position_vector(1, 0) + Scalar(110.0);
  const Scalar _tmp60 = _tmp58 + _tmp59;
  const Scalar _tmp61 = -_tmp27;
  const Scalar _tmp62 = _tmp30 + _tmp61;
  const Scalar _tmp63 = _tmp32 + _tmp62;
  const Scalar _tmp64 = std::pow(Scalar(std::pow(_tmp60, Scalar(2)) + std::pow(_tmp63, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp65 = _tmp63 * _tmp64;
  const Scalar _tmp66 = _tmp15 - _tmp7;
  const Scalar _tmp67 = _tmp60 * _tmp64;
  const Scalar _tmp68 = _tmp35 * _tmp65 - _tmp67;
  const Scalar _tmp69 = _tmp50 * _tmp68;
  const Scalar _tmp70 = _tmp16 * _tmp65 -
                        _tmp54 * (-_tmp36 * _tmp65 - _tmp49 * _tmp69 + _tmp66 * _tmp67) -
                        _tmp55 * _tmp69 - _tmp65 * _tmp66;
  const Scalar _tmp71 = Scalar(1.0) / (_tmp70);
  const Scalar _tmp72 =
      std::sqrt(Scalar(std::pow(_tmp26, Scalar(2)) + std::pow(_tmp33, Scalar(2))));
  const Scalar _tmp73 = Scalar(1.0) / (_tmp72);
  const Scalar _tmp74 = _tmp34 * _tmp72;
  const Scalar _tmp75 = _tmp74 * (-_tmp24 * _tmp33 * _tmp73 + _tmp26 * _tmp31 * _tmp73);
  const Scalar _tmp76 = _tmp38 * _tmp45 - _tmp41 * _tmp48 + _tmp45 * _tmp75;
  const Scalar _tmp77 = _tmp58 * _tmp65 - _tmp62 * _tmp67 + _tmp65 * _tmp75 - _tmp69 * _tmp76;
  const Scalar _tmp78 = _tmp71 * _tmp77;
  const Scalar _tmp79 = Scalar(1.0) / (_tmp77);
  const Scalar _tmp80 = _tmp70 * _tmp79;
  const Scalar _tmp81 = _tmp56 + _tmp80 * (_tmp51 * _tmp76 - _tmp56 * _tmp78 - _tmp75);
  const Scalar _tmp82 = _tmp68 * _tmp71;
  const Scalar _tmp83 = -_tmp35 - _tmp81 * _tmp82;
  const Scalar _tmp84 = _tmp45 * _tmp50;
  const Scalar _tmp85 = _tmp65 * _tmp71;
  const Scalar _tmp86 = _tmp40 + _tmp61;
  const Scalar _tmp87 = _tmp42 + _tmp86;
  const Scalar _tmp88 = _tmp37 + _tmp57;
  const Scalar _tmp89 = _tmp59 + _tmp88;
  const Scalar _tmp90 = std::pow(Scalar(std::pow(_tmp87, Scalar(2)) + std::pow(_tmp89, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp91 = _tmp87 * _tmp90;
  const Scalar _tmp92 = _tmp74 * fh1;
  const Scalar _tmp93 = Scalar(1.0) * _tmp50;
  const Scalar _tmp94 = _tmp49 * _tmp54 * _tmp93 - _tmp55 * _tmp93;
  const Scalar _tmp95 = _tmp80 * (-_tmp76 * _tmp93 - _tmp78 * _tmp94) + _tmp94;
  const Scalar _tmp96 = -_tmp82 * _tmp95 + Scalar(1.0);
  const Scalar _tmp97 = _tmp89 * _tmp90;
  const Scalar _tmp98 = Scalar(1.0) * _tmp79;
  const Scalar _tmp99 = fh1 * (_tmp86 * _tmp97 - _tmp88 * _tmp91);
  const Scalar _tmp100 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp101 = _tmp52 * _tmp54 + _tmp53;
  const Scalar _tmp102 = 0;
  const Scalar _tmp103 = _tmp102 * _tmp69;
  const Scalar _tmp104 = _tmp74 * (_tmp102 * _tmp65 - _tmp103 * _tmp45);
  const Scalar _tmp105 = _tmp98 * _tmp99;
  const Scalar _tmp106 = _tmp100 * _tmp102;
  const Scalar _tmp107 = _tmp50 * fh1;
  const Scalar _tmp108 = _tmp71 * fh1;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      -_tmp104 *
      std::exp(_tmp100 * _tmp104 + _tmp74 * _tmp99 * (-_tmp45 * _tmp69 * _tmp98 + _tmp65 * _tmp98) +
               _tmp91 * _tmp92 * (_tmp81 * _tmp85 + _tmp83 * _tmp84 + Scalar(1.0)) +
               _tmp92 * _tmp97 * (_tmp84 * _tmp96 + _tmp85 * _tmp95));
  _res(2, 0) = -_tmp103 * std::exp(_tmp105 * _tmp69 + _tmp106 * _tmp69 - _tmp107 * _tmp83 * _tmp91 -
                                   _tmp107 * _tmp96 * _tmp97);
  _res(3, 0) = _tmp102 *
               std::exp(-_tmp105 - _tmp106 - _tmp108 * _tmp81 * _tmp91 - _tmp108 * _tmp95 * _tmp97);

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
