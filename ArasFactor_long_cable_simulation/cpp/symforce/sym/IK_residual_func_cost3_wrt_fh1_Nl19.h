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
 * Symbolic function: IK_residual_func_cost3_wrt_fh1_Nl19
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost3WrtFh1Nl19(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 295

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (112)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp1 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp2 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp3 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = -Scalar(0.010999999999999999) * _tmp1 - Scalar(0.010999999999999999) * _tmp3;
  const Scalar _tmp5 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp6 = 2 * _tmp0;
  const Scalar _tmp7 = _tmp5 * _tmp6;
  const Scalar _tmp8 = -_DeltaRot[0] * _Rot_init[0] - _DeltaRot[1] * _Rot_init[1] -
                       _DeltaRot[2] * _Rot_init[2] + _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp9 = 2 * _tmp2 * _tmp8;
  const Scalar _tmp10 = Scalar(0.20999999999999999) * _tmp7 - Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp11 = 2 * _tmp5;
  const Scalar _tmp12 = _tmp11 * _tmp2;
  const Scalar _tmp13 = _tmp6 * _tmp8;
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp10 + _tmp14 + _tmp4;
  const Scalar _tmp16 = _tmp2 * _tmp6;
  const Scalar _tmp17 = _tmp11 * _tmp8;
  const Scalar _tmp18 = Scalar(0.20999999999999999) * _tmp16 + Scalar(0.20999999999999999) * _tmp17;
  const Scalar _tmp19 =
      -Scalar(0.010999999999999999) * _tmp12 + Scalar(0.010999999999999999) * _tmp13;
  const Scalar _tmp20 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp20 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp22 = _tmp19 + _tmp21;
  const Scalar _tmp23 = _tmp18 + _tmp22;
  const Scalar _tmp24 = position_vector(1, 0) + Scalar(-110.0);
  const Scalar _tmp25 = _tmp23 + _tmp24;
  const Scalar _tmp26 = Scalar(0.20999999999999999) * _tmp16 - Scalar(0.20999999999999999) * _tmp17;
  const Scalar _tmp27 =
      -Scalar(0.010999999999999999) * _tmp7 - Scalar(0.010999999999999999) * _tmp9;
  const Scalar _tmp28 = Scalar(0.20999999999999999) * _tmp20 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp29 = _tmp27 + _tmp28;
  const Scalar _tmp30 = _tmp26 + _tmp29;
  const Scalar _tmp31 = position_vector(0, 0) + Scalar(-125.0);
  const Scalar _tmp32 = _tmp30 + _tmp31;
  const Scalar _tmp33 = std::pow(Scalar(std::pow(_tmp25, Scalar(2)) + std::pow(_tmp32, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp34 = _tmp25 * _tmp33;
  const Scalar _tmp35 = -_tmp14 + _tmp4;
  const Scalar _tmp36 = -_tmp10 + _tmp35;
  const Scalar _tmp37 = -_tmp18;
  const Scalar _tmp38 = _tmp19 - _tmp21;
  const Scalar _tmp39 = _tmp37 + _tmp38;
  const Scalar _tmp40 = position_vector(1, 0) + Scalar(110.0);
  const Scalar _tmp41 = _tmp39 + _tmp40;
  const Scalar _tmp42 = -_tmp26;
  const Scalar _tmp43 = _tmp27 - _tmp28;
  const Scalar _tmp44 = _tmp42 + _tmp43;
  const Scalar _tmp45 = position_vector(0, 0) + Scalar(125.0);
  const Scalar _tmp46 = _tmp44 + _tmp45;
  const Scalar _tmp47 = Scalar(1.0) / (_tmp46);
  const Scalar _tmp48 = _tmp41 * _tmp47;
  const Scalar _tmp49 = _tmp36 * _tmp48;
  const Scalar _tmp50 = _tmp32 * _tmp33;
  const Scalar _tmp51 = _tmp15 * _tmp34 - _tmp49 * _tmp50;
  const Scalar _tmp52 = Scalar(1.0) * _tmp44;
  const Scalar _tmp53 = Scalar(1.0) * _tmp39;
  const Scalar _tmp54 = (-_tmp30 + _tmp52) / (_tmp23 - _tmp53);
  const Scalar _tmp55 = Scalar(1.0) / (-_tmp34 + _tmp48 * _tmp50);
  const Scalar _tmp56 = Scalar(1.0) * _tmp55;
  const Scalar _tmp57 = -_tmp15 * _tmp50 + _tmp36 * _tmp50;
  const Scalar _tmp58 = _tmp51 * _tmp54 * _tmp56 - _tmp56 * _tmp57;
  const Scalar _tmp59 =
      std::sqrt(Scalar(std::pow(_tmp41, Scalar(2)) + std::pow(_tmp46, Scalar(2))));
  const Scalar _tmp60 = Scalar(1.0) / (_tmp59);
  const Scalar _tmp61 = _tmp47 * _tmp59;
  const Scalar _tmp62 = _tmp61 * (-_tmp39 * _tmp46 * _tmp60 + _tmp41 * _tmp44 * _tmp60);
  const Scalar _tmp63 = _tmp23 * _tmp50 - _tmp30 * _tmp34 + _tmp50 * _tmp62;
  const Scalar _tmp64 = _tmp18 + _tmp38;
  const Scalar _tmp65 = _tmp40 + _tmp64;
  const Scalar _tmp66 = _tmp29 + _tmp42;
  const Scalar _tmp67 = _tmp31 + _tmp66;
  const Scalar _tmp68 = std::pow(Scalar(std::pow(_tmp65, Scalar(2)) + std::pow(_tmp67, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp69 = _tmp65 * _tmp68;
  const Scalar _tmp70 = _tmp67 * _tmp68;
  const Scalar _tmp71 = _tmp48 * _tmp70 - _tmp69;
  const Scalar _tmp72 = _tmp55 * _tmp71;
  const Scalar _tmp73 = _tmp62 * _tmp70 - _tmp63 * _tmp72 + _tmp64 * _tmp70 - _tmp66 * _tmp69;
  const Scalar _tmp74 = _tmp10 + _tmp35;
  const Scalar _tmp75 = _tmp36 * _tmp70;
  const Scalar _tmp76 = -_tmp54 * (-_tmp48 * _tmp75 - _tmp51 * _tmp72 + _tmp69 * _tmp74) -
                        _tmp57 * _tmp72 - _tmp70 * _tmp74 + _tmp75;
  const Scalar _tmp77 = Scalar(1.0) / (_tmp76);
  const Scalar _tmp78 = _tmp73 * _tmp77;
  const Scalar _tmp79 = Scalar(1.0) / (_tmp73);
  const Scalar _tmp80 = _tmp76 * _tmp79;
  const Scalar _tmp81 = _tmp58 + _tmp80 * (-_tmp56 * _tmp63 - _tmp58 * _tmp78);
  const Scalar _tmp82 = _tmp70 * _tmp77;
  const Scalar _tmp83 = _tmp71 * _tmp77;
  const Scalar _tmp84 = _tmp55 * (-_tmp81 * _tmp83 + Scalar(1.0));
  const Scalar _tmp85 = _tmp26 + _tmp43;
  const Scalar _tmp86 = _tmp45 + _tmp85;
  const Scalar _tmp87 = _tmp22 + _tmp37;
  const Scalar _tmp88 = _tmp24 + _tmp87;
  const Scalar _tmp89 = std::pow(Scalar(std::pow(_tmp86, Scalar(2)) + std::pow(_tmp88, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp90 = _tmp88 * _tmp89;
  const Scalar _tmp91 = _tmp61 * _tmp90 * (_tmp50 * _tmp84 + _tmp81 * _tmp82);
  const Scalar _tmp92 = _tmp86 * _tmp89;
  const Scalar _tmp93 = _tmp85 * _tmp90 - _tmp87 * _tmp92;
  const Scalar _tmp94 = Scalar(1.0) * _tmp79;
  const Scalar _tmp95 = _tmp50 * _tmp72;
  const Scalar _tmp96 = _tmp61 * _tmp93 * (_tmp70 * _tmp94 - _tmp94 * _tmp95);
  const Scalar _tmp97 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp98 = _tmp52 + _tmp53 * _tmp54;
  const Scalar _tmp99 = 0;
  const Scalar _tmp100 = _tmp48 * _tmp55;
  const Scalar _tmp101 = _tmp100 * _tmp57 - _tmp36 - _tmp54 * (_tmp100 * _tmp51 + _tmp49);
  const Scalar _tmp102 = _tmp101 + _tmp80 * (_tmp100 * _tmp63 - _tmp101 * _tmp78 - _tmp62);
  const Scalar _tmp103 = _tmp55 * (-_tmp102 * _tmp83 - _tmp48);
  const Scalar _tmp104 = _tmp61 * _tmp92 * (_tmp102 * _tmp82 + _tmp103 * _tmp50 + Scalar(1.0));
  const Scalar _tmp105 = _tmp103 * _tmp92;
  const Scalar _tmp106 = _tmp93 * _tmp94;
  const Scalar _tmp107 = _tmp106 * fh1;
  const Scalar _tmp108 = _tmp84 * _tmp90;
  const Scalar _tmp109 = _tmp97 * _tmp99;
  const Scalar _tmp110 = _tmp77 * _tmp81 * _tmp90;
  const Scalar _tmp111 = _tmp102 * _tmp77 * _tmp92;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = -std::exp(-fh1);
  _res(1, 0) = -(-_tmp104 - _tmp91 - _tmp96) *
               std::exp(_tmp104 * fh1 + _tmp61 * _tmp97 * (_tmp70 * _tmp99 - _tmp95 * _tmp99) +
                        _tmp91 * fh1 + _tmp96 * fh1);
  _res(2, 0) = -(_tmp105 - _tmp106 * _tmp72 + _tmp108) *
               std::exp(-_tmp105 * fh1 + _tmp107 * _tmp72 - _tmp108 * fh1 + _tmp109 * _tmp72);
  _res(3, 0) =
      -(_tmp106 + _tmp110 + _tmp111) * std::exp(-_tmp107 - _tmp109 - _tmp110 * fh1 - _tmp111 * fh1);

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
