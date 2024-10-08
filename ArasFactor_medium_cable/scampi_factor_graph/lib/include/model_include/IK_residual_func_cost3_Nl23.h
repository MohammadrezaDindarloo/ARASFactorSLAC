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
 * Symbolic function: IK_residual_func_cost3_Nl23
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost3Nl23(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 285

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (100)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp2 = 2 * _tmp0 * _tmp1;
  const Scalar _tmp3 = -_DeltaRot[0] * _Rot_init[0] - _DeltaRot[1] * _Rot_init[1] -
                       _DeltaRot[2] * _Rot_init[2] + _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp4 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp5 = 2 * _tmp4;
  const Scalar _tmp6 = _tmp3 * _tmp5;
  const Scalar _tmp7 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp8 = _tmp0 * _tmp5;
  const Scalar _tmp9 = 2 * _tmp3;
  const Scalar _tmp10 = _tmp1 * _tmp9;
  const Scalar _tmp11 =
      Scalar(0.010999999999999999) * _tmp10 - Scalar(0.010999999999999999) * _tmp8;
  const Scalar _tmp12 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp13 = 1 - 2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp11 + _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp7;
  const Scalar _tmp17 = _tmp16 + position_vector(1, 0) + Scalar(-4.7744369927459998);
  const Scalar _tmp18 = Scalar(0.20999999999999999) * _tmp2 - Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp19 = _tmp1 * _tmp5;
  const Scalar _tmp20 = _tmp0 * _tmp9;
  const Scalar _tmp21 =
      -Scalar(0.010999999999999999) * _tmp19 - Scalar(0.010999999999999999) * _tmp20;
  const Scalar _tmp22 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp23 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp22;
  const Scalar _tmp24 = _tmp21 + _tmp23;
  const Scalar _tmp25 = _tmp18 + _tmp24;
  const Scalar _tmp26 = _tmp25 + position_vector(0, 0) + Scalar(-2.7171519410699099);
  const Scalar _tmp27 = Scalar(1.0) / (_tmp26);
  const Scalar _tmp28 = _tmp17 * _tmp27;
  const Scalar _tmp29 = _tmp11 - _tmp14;
  const Scalar _tmp30 = _tmp29 + _tmp7;
  const Scalar _tmp31 = _tmp30 + position_vector(1, 0) + Scalar(8.3885017487099702);
  const Scalar _tmp32 = -_tmp18;
  const Scalar _tmp33 = _tmp24 + _tmp32;
  const Scalar _tmp34 = _tmp33 + position_vector(0, 0) + Scalar(-2.5193355532036801);
  const Scalar _tmp35 = std::pow(Scalar(std::pow(_tmp31, Scalar(2)) + std::pow(_tmp34, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp36 = _tmp34 * _tmp35;
  const Scalar _tmp37 = _tmp31 * _tmp35;
  const Scalar _tmp38 = Scalar(1.0) / (_tmp28 * _tmp36 - _tmp37);
  const Scalar _tmp39 =
      std::sqrt(Scalar(std::pow(_tmp17, Scalar(2)) + std::pow(_tmp26, Scalar(2))));
  const Scalar _tmp40 = Scalar(1.0) / (_tmp39);
  const Scalar _tmp41 = _tmp27 * _tmp39;
  const Scalar _tmp42 = _tmp41 * (-_tmp16 * _tmp26 * _tmp40 + _tmp17 * _tmp25 * _tmp40);
  const Scalar _tmp43 = _tmp38 * (_tmp30 * _tmp36 - _tmp33 * _tmp37 + _tmp36 * _tmp42);
  const Scalar _tmp44 = -Scalar(0.010999999999999999) * _tmp12 -
                        Scalar(0.010999999999999999) * _tmp22 + Scalar(-0.010999999999999999);
  const Scalar _tmp45 = Scalar(0.20999999999999999) * _tmp19 - Scalar(0.20999999999999999) * _tmp20;
  const Scalar _tmp46 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp47 = _tmp44 + _tmp45 + _tmp46;
  const Scalar _tmp48 = _tmp44 - _tmp46;
  const Scalar _tmp49 = _tmp45 + _tmp48;
  const Scalar _tmp50 = _tmp28 * _tmp47;
  const Scalar _tmp51 = _tmp38 * (-_tmp36 * _tmp50 + _tmp37 * _tmp49);
  const Scalar _tmp52 = Scalar(1.0) * _tmp16;
  const Scalar _tmp53 = Scalar(1.0) * _tmp25;
  const Scalar _tmp54 = (-_tmp33 + _tmp53) / (_tmp30 - _tmp52);
  const Scalar _tmp55 = _tmp38 * (_tmp36 * _tmp47 - _tmp36 * _tmp49);
  const Scalar _tmp56 = _tmp28 * _tmp55 - _tmp47 - _tmp54 * (_tmp28 * _tmp51 + _tmp50);
  const Scalar _tmp57 = -_tmp45 + _tmp48;
  const Scalar _tmp58 = -_tmp7;
  const Scalar _tmp59 = _tmp29 + _tmp58;
  const Scalar _tmp60 = _tmp59 + position_vector(1, 0) + Scalar(8.3196563720703107);
  const Scalar _tmp61 = _tmp21 - _tmp23;
  const Scalar _tmp62 = _tmp32 + _tmp61;
  const Scalar _tmp63 = _tmp62 + position_vector(0, 0) + Scalar(1.9874742031097401);
  const Scalar _tmp64 = std::pow(Scalar(std::pow(_tmp60, Scalar(2)) + std::pow(_tmp63, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp65 = _tmp60 * _tmp64;
  const Scalar _tmp66 = _tmp63 * _tmp64;
  const Scalar _tmp67 = _tmp28 * _tmp66 - _tmp65;
  const Scalar _tmp68 = _tmp47 * _tmp66 -
                        _tmp54 * (-_tmp50 * _tmp66 - _tmp51 * _tmp67 + _tmp57 * _tmp65) -
                        _tmp55 * _tmp67 - _tmp57 * _tmp66;
  const Scalar _tmp69 = Scalar(1.0) / (_tmp68);
  const Scalar _tmp70 = _tmp42 * _tmp66 - _tmp43 * _tmp67 + _tmp59 * _tmp66 - _tmp62 * _tmp65;
  const Scalar _tmp71 = _tmp69 * _tmp70;
  const Scalar _tmp72 = Scalar(1.0) / (_tmp70);
  const Scalar _tmp73 = _tmp68 * _tmp72;
  const Scalar _tmp74 = _tmp56 + _tmp73 * (_tmp28 * _tmp43 - _tmp42 - _tmp56 * _tmp71);
  const Scalar _tmp75 = _tmp67 * _tmp69;
  const Scalar _tmp76 = -_tmp28 - _tmp74 * _tmp75;
  const Scalar _tmp77 = _tmp36 * _tmp38;
  const Scalar _tmp78 = _tmp66 * _tmp69;
  const Scalar _tmp79 = _tmp18 + _tmp61;
  const Scalar _tmp80 = _tmp79 + position_vector(0, 0) + Scalar(1.7965602546229);
  const Scalar _tmp81 = _tmp15 + _tmp58;
  const Scalar _tmp82 = _tmp81 + position_vector(1, 0) + Scalar(-4.83288938413423);
  const Scalar _tmp83 = std::pow(Scalar(std::pow(_tmp80, Scalar(2)) + std::pow(_tmp82, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp84 = _tmp80 * _tmp83;
  const Scalar _tmp85 = _tmp41 * fh1;
  const Scalar _tmp86 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp87 = _tmp52 * _tmp54 + _tmp53;
  const Scalar _tmp88 = 0;
  const Scalar _tmp89 = _tmp69 * _tmp88;
  const Scalar _tmp90 = _tmp75 * _tmp88;
  const Scalar _tmp91 = Scalar(1.0) * _tmp51 * _tmp54 - Scalar(1.0) * _tmp55;
  const Scalar _tmp92 = _tmp73 * (-Scalar(1.0) * _tmp43 - _tmp71 * _tmp91) + _tmp91;
  const Scalar _tmp93 = -_tmp75 * _tmp92 + Scalar(1.0);
  const Scalar _tmp94 = _tmp82 * _tmp83;
  const Scalar _tmp95 = Scalar(1.0) * _tmp72;
  const Scalar _tmp96 = fh1 * (_tmp79 * _tmp94 - _tmp81 * _tmp84);
  const Scalar _tmp97 = _tmp95 * _tmp96;
  const Scalar _tmp98 = _tmp38 * fh1;
  const Scalar _tmp99 = _tmp69 * fh1;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = std::exp(-fh1);
  _res(1, 0) = std::exp(_tmp41 * _tmp86 * (_tmp66 * _tmp89 - _tmp77 * _tmp90) +
                        _tmp41 * _tmp96 * (_tmp66 * _tmp95 - _tmp67 * _tmp77 * _tmp95) +
                        _tmp84 * _tmp85 * (_tmp74 * _tmp78 + _tmp76 * _tmp77 + Scalar(1.0)) +
                        _tmp85 * _tmp94 * (_tmp77 * _tmp93 + _tmp78 * _tmp92));
  _res(2, 0) = std::exp(_tmp38 * _tmp67 * _tmp97 + _tmp38 * _tmp86 * _tmp90 -
                        _tmp76 * _tmp84 * _tmp98 - _tmp93 * _tmp94 * _tmp98);
  _res(3, 0) =
      std::exp(-_tmp74 * _tmp84 * _tmp99 - _tmp86 * _tmp89 - _tmp92 * _tmp94 * _tmp99 - _tmp97);

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
