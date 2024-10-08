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
 * Symbolic function: IK_residual_func_cost3_wrt_fv1_Nl7
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost3WrtFv1Nl7(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 284

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (108)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp2 = 2 * _tmp1;
  const Scalar _tmp3 = _tmp0 * _tmp2;
  const Scalar _tmp4 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp5 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp6 = _tmp4 * _tmp5;
  const Scalar _tmp7 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp8 = -_tmp7;
  const Scalar _tmp9 = _tmp2 * _tmp4;
  const Scalar _tmp10 = _tmp0 * _tmp5;
  const Scalar _tmp11 =
      Scalar(0.010999999999999999) * _tmp10 - Scalar(0.010999999999999999) * _tmp9;
  const Scalar _tmp12 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp13 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 +
                        Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999);
  const Scalar _tmp15 = _tmp11 - _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp8;
  const Scalar _tmp17 = position_vector(1, 0) + Scalar(110.0);
  const Scalar _tmp18 = _tmp16 + _tmp17;
  const Scalar _tmp19 = Scalar(0.20999999999999999) * _tmp3 - Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp20 = -_tmp19;
  const Scalar _tmp21 = 2 * _tmp0 * _tmp4;
  const Scalar _tmp22 = _tmp1 * _tmp5;
  const Scalar _tmp23 =
      -Scalar(0.010999999999999999) * _tmp21 - Scalar(0.010999999999999999) * _tmp22;
  const Scalar _tmp24 = 1 - 2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp25 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp26 = _tmp23 - _tmp25;
  const Scalar _tmp27 = _tmp20 + _tmp26;
  const Scalar _tmp28 = position_vector(0, 0) + Scalar(125.0);
  const Scalar _tmp29 = _tmp27 + _tmp28;
  const Scalar _tmp30 = Scalar(1.0) / (_tmp29);
  const Scalar _tmp31 = _tmp18 * _tmp30;
  const Scalar _tmp32 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp33 =
      -Scalar(0.010999999999999999) * _tmp12 - Scalar(0.010999999999999999) * _tmp24;
  const Scalar _tmp34 = Scalar(0.20999999999999999) * _tmp21 - Scalar(0.20999999999999999) * _tmp22;
  const Scalar _tmp35 = _tmp33 - _tmp34;
  const Scalar _tmp36 = -_tmp32 + _tmp35;
  const Scalar _tmp37 = _tmp19 + _tmp26;
  const Scalar _tmp38 = _tmp28 + _tmp37;
  const Scalar _tmp39 = _tmp11 + _tmp14;
  const Scalar _tmp40 = _tmp39 + _tmp8;
  const Scalar _tmp41 = position_vector(1, 0) + Scalar(-110.0);
  const Scalar _tmp42 = _tmp40 + _tmp41;
  const Scalar _tmp43 = std::pow(Scalar(std::pow(_tmp38, Scalar(2)) + std::pow(_tmp42, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp44 = _tmp38 * _tmp43;
  const Scalar _tmp45 = _tmp36 * _tmp44;
  const Scalar _tmp46 = _tmp32 + _tmp35;
  const Scalar _tmp47 = _tmp42 * _tmp43;
  const Scalar _tmp48 = -_tmp31 * _tmp45 + _tmp46 * _tmp47;
  const Scalar _tmp49 = Scalar(1.0) * _tmp16;
  const Scalar _tmp50 = Scalar(1.0) * _tmp27;
  const Scalar _tmp51 = (-_tmp37 + _tmp50) / (_tmp40 - _tmp49);
  const Scalar _tmp52 = Scalar(1.0) / (_tmp31 * _tmp44 - _tmp47);
  const Scalar _tmp53 = Scalar(1.0) * _tmp52;
  const Scalar _tmp54 = -_tmp44 * _tmp46 + _tmp45;
  const Scalar _tmp55 = _tmp48 * _tmp51 * _tmp53 - _tmp53 * _tmp54;
  const Scalar _tmp56 = _tmp39 + _tmp7;
  const Scalar _tmp57 = _tmp41 + _tmp56;
  const Scalar _tmp58 = _tmp23 + _tmp25;
  const Scalar _tmp59 = _tmp19 + _tmp58;
  const Scalar _tmp60 = position_vector(0, 0) + Scalar(-125.0);
  const Scalar _tmp61 = _tmp59 + _tmp60;
  const Scalar _tmp62 = std::pow(Scalar(std::pow(_tmp57, Scalar(2)) + std::pow(_tmp61, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp63 = _tmp61 * _tmp62;
  const Scalar _tmp64 = _tmp36 * _tmp63;
  const Scalar _tmp65 = _tmp32 + _tmp33 + _tmp34;
  const Scalar _tmp66 = _tmp57 * _tmp62;
  const Scalar _tmp67 = _tmp31 * _tmp63 - _tmp66;
  const Scalar _tmp68 = _tmp52 * _tmp67;
  const Scalar _tmp69 = -_tmp51 * (-_tmp31 * _tmp64 - _tmp48 * _tmp68 + _tmp65 * _tmp66) -
                        _tmp54 * _tmp68 - _tmp63 * _tmp65 + _tmp64;
  const Scalar _tmp70 = Scalar(1.0) / (_tmp69);
  const Scalar _tmp71 =
      std::sqrt(Scalar(std::pow(_tmp18, Scalar(2)) + std::pow(_tmp29, Scalar(2))));
  const Scalar _tmp72 = Scalar(1.0) / (_tmp71);
  const Scalar _tmp73 = _tmp30 * _tmp71;
  const Scalar _tmp74 = _tmp73 * (-_tmp16 * _tmp29 * _tmp72 + _tmp18 * _tmp27 * _tmp72);
  const Scalar _tmp75 = -_tmp37 * _tmp47 + _tmp40 * _tmp44 + _tmp44 * _tmp74;
  const Scalar _tmp76 = _tmp56 * _tmp63 - _tmp59 * _tmp66 + _tmp63 * _tmp74 - _tmp68 * _tmp75;
  const Scalar _tmp77 = _tmp70 * _tmp76;
  const Scalar _tmp78 = Scalar(1.0) / (_tmp76);
  const Scalar _tmp79 = _tmp69 * _tmp78;
  const Scalar _tmp80 = _tmp55 + _tmp79 * (-_tmp53 * _tmp75 - _tmp55 * _tmp77);
  const Scalar _tmp81 = _tmp63 * _tmp70;
  const Scalar _tmp82 = _tmp67 * _tmp70;
  const Scalar _tmp83 = -_tmp80 * _tmp82 + Scalar(1.0);
  const Scalar _tmp84 = _tmp44 * _tmp52;
  const Scalar _tmp85 = _tmp15 + _tmp7;
  const Scalar _tmp86 = _tmp17 + _tmp85;
  const Scalar _tmp87 = _tmp20 + _tmp58;
  const Scalar _tmp88 = _tmp60 + _tmp87;
  const Scalar _tmp89 = std::pow(Scalar(std::pow(_tmp86, Scalar(2)) + std::pow(_tmp88, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp90 = _tmp86 * _tmp89;
  const Scalar _tmp91 = _tmp73 * fh1;
  const Scalar _tmp92 = Scalar(1.0) * _tmp78;
  const Scalar _tmp93 = _tmp88 * _tmp89;
  const Scalar _tmp94 = fh1 * (-_tmp85 * _tmp93 + _tmp87 * _tmp90);
  const Scalar _tmp95 = _tmp31 * _tmp52;
  const Scalar _tmp96 = -_tmp36 - _tmp51 * (_tmp31 * _tmp36 + _tmp48 * _tmp95) + _tmp54 * _tmp95;
  const Scalar _tmp97 = _tmp79 * (-_tmp74 + _tmp75 * _tmp95 - _tmp77 * _tmp96) + _tmp96;
  const Scalar _tmp98 = -_tmp31 - _tmp82 * _tmp97;
  const Scalar _tmp99 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp100 = _tmp49 * _tmp51 + _tmp50;
  const Scalar _tmp101 = 0;
  const Scalar _tmp102 = _tmp101 * _tmp68;
  const Scalar _tmp103 = _tmp73 * (_tmp101 * _tmp63 - _tmp102 * _tmp44);
  const Scalar _tmp104 = _tmp52 * fh1;
  const Scalar _tmp105 = _tmp101 * _tmp99;
  const Scalar _tmp106 = _tmp92 * _tmp94;
  const Scalar _tmp107 = _tmp70 * fh1;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      -_tmp103 *
      std::exp(_tmp103 * _tmp99 + _tmp73 * _tmp94 * (-_tmp44 * _tmp68 * _tmp92 + _tmp63 * _tmp92) +
               _tmp90 * _tmp91 * (_tmp80 * _tmp81 + _tmp83 * _tmp84) +
               _tmp91 * _tmp93 * (_tmp81 * _tmp97 + _tmp84 * _tmp98 + Scalar(1.0)));
  _res(2, 0) = -_tmp102 * std::exp(-_tmp104 * _tmp83 * _tmp90 - _tmp104 * _tmp93 * _tmp98 +
                                   _tmp105 * _tmp68 + _tmp106 * _tmp68);
  _res(3, 0) = _tmp101 *
               std::exp(-_tmp105 - _tmp106 - _tmp107 * _tmp80 * _tmp90 - _tmp107 * _tmp93 * _tmp97);

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
