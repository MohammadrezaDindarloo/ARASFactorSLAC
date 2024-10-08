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
 * Symbolic function: IK_residual_func_cost3_Nl18
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     p_init0: Scalar
 *     p_init1: Scalar
 *     p_init2: Scalar
 *     rot_init_x: Scalar
 *     rot_init_y: Scalar
 *     rot_init_z: Scalar
 *     rot_init_w: Scalar
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix41
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost3Nl18(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot, const Scalar p_init0,
    const Scalar p_init1, const Scalar p_init2, const Scalar rot_init_x, const Scalar rot_init_y,
    const Scalar rot_init_z, const Scalar rot_init_w, const Scalar epsilon) {
  // Total ops: 288

  // Unused inputs
  (void)p_init2;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();

  // Intermediate terms (100)
  const Scalar _tmp0 = _DeltaRot[0] * rot_init_z + _DeltaRot[1] * rot_init_w -
                       _DeltaRot[2] * rot_init_x + _DeltaRot[3] * rot_init_y;
  const Scalar _tmp1 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp2 = -_DeltaRot[0] * rot_init_y + _DeltaRot[1] * rot_init_x +
                       _DeltaRot[2] * rot_init_w + _DeltaRot[3] * rot_init_z;
  const Scalar _tmp3 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp5 = _DeltaRot[0] * rot_init_w - _DeltaRot[1] * rot_init_z +
                       _DeltaRot[2] * rot_init_y + _DeltaRot[3] * rot_init_x;
  const Scalar _tmp6 = 2 * _tmp2 * _tmp5;
  const Scalar _tmp7 = -2 * _DeltaRot[0] * rot_init_x - 2 * _DeltaRot[1] * rot_init_y -
                       2 * _DeltaRot[2] * rot_init_z + 2 * _DeltaRot[3] * rot_init_w;
  const Scalar _tmp8 = _tmp0 * _tmp7;
  const Scalar _tmp9 = -Scalar(0.010999999999999999) * _tmp6 - Scalar(0.010999999999999999) * _tmp8;
  const Scalar _tmp10 = 2 * _tmp0;
  const Scalar _tmp11 = _tmp10 * _tmp5;
  const Scalar _tmp12 = _tmp2 * _tmp7;
  const Scalar _tmp13 = Scalar(0.20999999999999999) * _tmp11 - Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp14 = _tmp13 + _tmp9;
  const Scalar _tmp15 = _tmp14 + _tmp4;
  const Scalar _tmp16 = _tmp15 + p_init0 + Scalar(-2.71799795);
  const Scalar _tmp17 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp18 = Scalar(0.20999999999999999) * _tmp17 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp19 = _tmp10 * _tmp2;
  const Scalar _tmp20 = _tmp5 * _tmp7;
  const Scalar _tmp21 =
      -Scalar(0.010999999999999999) * _tmp19 + Scalar(0.010999999999999999) * _tmp20;
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp23 = _tmp21 + _tmp22;
  const Scalar _tmp24 = _tmp18 + _tmp23;
  const Scalar _tmp25 = _tmp24 + p_init1 + Scalar(-4.7752063900000001);
  const Scalar _tmp26 = std::pow(Scalar(std::pow(_tmp16, Scalar(2)) + std::pow(_tmp25, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp27 = _tmp16 * _tmp26;
  const Scalar _tmp28 = _tmp25 * _tmp26;
  const Scalar _tmp29 = -_tmp18;
  const Scalar _tmp30 = _tmp21 - _tmp22;
  const Scalar _tmp31 = _tmp29 + _tmp30;
  const Scalar _tmp32 = _tmp31 + p_init1 + Scalar(8.3196563700000006);
  const Scalar _tmp33 = -_tmp4;
  const Scalar _tmp34 = -_tmp13 + _tmp9;
  const Scalar _tmp35 = _tmp33 + _tmp34;
  const Scalar _tmp36 = _tmp35 + p_init0 + Scalar(1.9874742000000001);
  const Scalar _tmp37 =
      std::sqrt(Scalar(std::pow(_tmp32, Scalar(2)) + std::pow(_tmp36, Scalar(2))));
  const Scalar _tmp38 = Scalar(1.0) / (_tmp37);
  const Scalar _tmp39 = Scalar(1.0) / (_tmp36);
  const Scalar _tmp40 = _tmp37 * _tmp39;
  const Scalar _tmp41 = _tmp40 * (-_tmp31 * _tmp36 * _tmp38 + _tmp32 * _tmp35 * _tmp38);
  const Scalar _tmp42 = _tmp32 * _tmp39;
  const Scalar _tmp43 = _tmp27 * _tmp42 - _tmp28;
  const Scalar _tmp44 = _tmp34 + _tmp4;
  const Scalar _tmp45 = _tmp44 + p_init0 + Scalar(-2.5202214700000001);
  const Scalar _tmp46 = _tmp23 + _tmp29;
  const Scalar _tmp47 = _tmp46 + p_init1 + Scalar(8.3888750099999996);
  const Scalar _tmp48 = std::pow(Scalar(std::pow(_tmp45, Scalar(2)) + std::pow(_tmp47, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp49 = _tmp45 * _tmp48;
  const Scalar _tmp50 = _tmp47 * _tmp48;
  const Scalar _tmp51 = Scalar(1.0) / (_tmp42 * _tmp49 - _tmp50);
  const Scalar _tmp52 = _tmp51 * (_tmp41 * _tmp49 - _tmp44 * _tmp50 + _tmp46 * _tmp49);
  const Scalar _tmp53 = -_tmp15 * _tmp28 + _tmp24 * _tmp27 + _tmp27 * _tmp41 - _tmp43 * _tmp52;
  const Scalar _tmp54 = Scalar(1.0) / (_tmp53);
  const Scalar _tmp55 = Scalar(1.0) * _tmp54;
  const Scalar _tmp56 = _tmp49 * _tmp51;
  const Scalar _tmp57 = _tmp14 + _tmp33;
  const Scalar _tmp58 = _tmp18 + _tmp30;
  const Scalar _tmp59 = _tmp58 + p_init1 + Scalar(-4.8333311099999996);
  const Scalar _tmp60 = _tmp57 + p_init0 + Scalar(1.79662371);
  const Scalar _tmp61 = std::pow(Scalar(std::pow(_tmp59, Scalar(2)) + std::pow(_tmp60, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp62 = _tmp59 * _tmp61;
  const Scalar _tmp63 = _tmp60 * _tmp61;
  const Scalar _tmp64 = fh1 * (_tmp57 * _tmp62 - _tmp58 * _tmp63);
  const Scalar _tmp65 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp66 = Scalar(1.0) * _tmp35;
  const Scalar _tmp67 = Scalar(1.0) * _tmp31;
  const Scalar _tmp68 = (-_tmp44 + _tmp66) / (_tmp46 - _tmp67);
  const Scalar _tmp69 = _tmp66 + _tmp67 * _tmp68;
  const Scalar _tmp70 = 0;
  const Scalar _tmp71 = Scalar(0.20999999999999999) * _tmp6 - Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp72 = -Scalar(0.010999999999999999) * _tmp1 -
                        Scalar(0.010999999999999999) * _tmp17 + Scalar(-0.010999999999999999);
  const Scalar _tmp73 = Scalar(0.20999999999999999) * _tmp19 + Scalar(0.20999999999999999) * _tmp20;
  const Scalar _tmp74 = _tmp72 - _tmp73;
  const Scalar _tmp75 = _tmp71 + _tmp74;
  const Scalar _tmp76 = -_tmp71 + _tmp74;
  const Scalar _tmp77 = _tmp49 * _tmp76;
  const Scalar _tmp78 = _tmp51 * (-_tmp49 * _tmp75 + _tmp77);
  const Scalar _tmp79 = _tmp27 * _tmp76;
  const Scalar _tmp80 = _tmp71 + _tmp72 + _tmp73;
  const Scalar _tmp81 = _tmp51 * (-_tmp42 * _tmp77 + _tmp50 * _tmp75);
  const Scalar _tmp82 = -_tmp27 * _tmp80 - _tmp43 * _tmp78 -
                        _tmp68 * (_tmp28 * _tmp80 - _tmp42 * _tmp79 - _tmp43 * _tmp81) + _tmp79;
  const Scalar _tmp83 = Scalar(1.0) / (_tmp82);
  const Scalar _tmp84 = _tmp43 * _tmp83;
  const Scalar _tmp85 = _tmp70 * _tmp84;
  const Scalar _tmp86 = _tmp70 * _tmp83;
  const Scalar _tmp87 = _tmp42 * _tmp78 - _tmp68 * (_tmp42 * _tmp76 + _tmp42 * _tmp81) - _tmp76;
  const Scalar _tmp88 = _tmp53 * _tmp83;
  const Scalar _tmp89 = _tmp54 * _tmp82;
  const Scalar _tmp90 = _tmp87 + _tmp89 * (-_tmp41 + _tmp42 * _tmp52 - _tmp87 * _tmp88);
  const Scalar _tmp91 = -_tmp42 - _tmp84 * _tmp90;
  const Scalar _tmp92 = _tmp27 * _tmp83;
  const Scalar _tmp93 = _tmp40 * fh1;
  const Scalar _tmp94 = Scalar(1.0) * _tmp68 * _tmp81 - Scalar(1.0) * _tmp78;
  const Scalar _tmp95 = _tmp89 * (-Scalar(1.0) * _tmp52 - _tmp88 * _tmp94) + _tmp94;
  const Scalar _tmp96 = -_tmp84 * _tmp95 + Scalar(1.0);
  const Scalar _tmp97 = _tmp55 * _tmp64;
  const Scalar _tmp98 = _tmp51 * fh1;
  const Scalar _tmp99 = _tmp83 * fh1;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = std::exp(-fh1);
  _res(1, 0) = std::exp(_tmp40 * _tmp64 * (_tmp27 * _tmp55 - _tmp43 * _tmp55 * _tmp56) +
                        _tmp40 * _tmp65 * (_tmp27 * _tmp86 - _tmp56 * _tmp85) +
                        _tmp62 * _tmp93 * (_tmp56 * _tmp96 + _tmp92 * _tmp95) +
                        _tmp63 * _tmp93 * (_tmp56 * _tmp91 + _tmp90 * _tmp92 + Scalar(1.0)));
  _res(2, 0) = std::exp(_tmp43 * _tmp51 * _tmp97 + _tmp51 * _tmp65 * _tmp85 -
                        _tmp62 * _tmp96 * _tmp98 - _tmp63 * _tmp91 * _tmp98);
  _res(3, 0) =
      std::exp(-_tmp62 * _tmp95 * _tmp99 - _tmp63 * _tmp90 * _tmp99 - _tmp65 * _tmp86 - _tmp97);

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
