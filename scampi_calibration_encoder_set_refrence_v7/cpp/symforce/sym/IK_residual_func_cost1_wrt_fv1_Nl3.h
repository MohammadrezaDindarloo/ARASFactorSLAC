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
 * Symbolic function: IK_residual_func_cost1_wrt_fv1_Nl3
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     p_a: Matrix31
 *     p_b: Matrix31
 *     p_c: Matrix31
 *     p_d: Matrix31
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1WrtFv1Nl3(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const Scalar p_init0, const Scalar p_init1, const Scalar p_init2, const Scalar rot_init_x,
    const Scalar rot_init_y, const Scalar rot_init_z, const Scalar rot_init_w,
    const Scalar epsilon) {
  // Total ops: 583

  // Unused inputs
  (void)p_init2;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();

  // Intermediate terms (192)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp2 = Scalar(1.0) * _tmp0 /
                       std::sqrt(Scalar(1 + std::pow(fv1, Scalar(2)) / std::pow(fh1, Scalar(2))));
  const Scalar _tmp3 = Scalar(9.6622558468725703) * fh1;
  const Scalar _tmp4 = _DeltaRot[0] * rot_init_w - _DeltaRot[1] * rot_init_z +
                       _DeltaRot[2] * rot_init_y + _DeltaRot[3] * rot_init_x;
  const Scalar _tmp5 = _DeltaRot[0] * rot_init_z + _DeltaRot[1] * rot_init_w -
                       _DeltaRot[2] * rot_init_x + _DeltaRot[3] * rot_init_y;
  const Scalar _tmp6 = 2 * _tmp5;
  const Scalar _tmp7 = _tmp4 * _tmp6;
  const Scalar _tmp8 = -_DeltaRot[0] * rot_init_y + _DeltaRot[1] * rot_init_x +
                       _DeltaRot[2] * rot_init_w + _DeltaRot[3] * rot_init_z;
  const Scalar _tmp9 = -2 * _DeltaRot[0] * rot_init_x - 2 * _DeltaRot[1] * rot_init_y -
                       2 * _DeltaRot[2] * rot_init_z + 2 * _DeltaRot[3] * rot_init_w;
  const Scalar _tmp10 = _tmp8 * _tmp9;
  const Scalar _tmp11 = -Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp12 = -_tmp11;
  const Scalar _tmp13 = 2 * _tmp4 * _tmp8;
  const Scalar _tmp14 = _tmp5 * _tmp9;
  const Scalar _tmp15 = _tmp13 + _tmp14;
  const Scalar _tmp16 = -Scalar(0.010999999999999999) * _tmp15;
  const Scalar _tmp17 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp18 = 1 - 2 * std::pow(_tmp8, Scalar(2));
  const Scalar _tmp19 = Scalar(0.20999999999999999) * _tmp17 + Scalar(0.20999999999999999) * _tmp18;
  const Scalar _tmp20 = _tmp16 - _tmp19;
  const Scalar _tmp21 = _tmp12 + _tmp20;
  const Scalar _tmp22 = _tmp21 + p_init0;
  const Scalar _tmp23 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp24 = Scalar(0.20999999999999999) * _tmp18 + Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp25 = -_tmp24;
  const Scalar _tmp26 = _tmp6 * _tmp8;
  const Scalar _tmp27 = _tmp4 * _tmp9;
  const Scalar _tmp28 = _tmp26 - _tmp27;
  const Scalar _tmp29 = Scalar(0.010999999999999999) * _tmp28;
  const Scalar _tmp30 = -_tmp29;
  const Scalar _tmp31 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp32 = _tmp30 - _tmp31;
  const Scalar _tmp33 = _tmp25 + _tmp32;
  const Scalar _tmp34 = _tmp33 + p_init1;
  const Scalar _tmp35 = _tmp30 + _tmp31;
  const Scalar _tmp36 = _tmp25 + _tmp35;
  const Scalar _tmp37 = _tmp36 + p_init1;
  const Scalar _tmp38 = _tmp37 - p_b(1, 0);
  const Scalar _tmp39 = _tmp16 + _tmp19;
  const Scalar _tmp40 = _tmp12 + _tmp39;
  const Scalar _tmp41 = _tmp40 + p_init0;
  const Scalar _tmp42 = _tmp41 - p_b(0, 0);
  const Scalar _tmp43 = std::pow(Scalar(std::pow(_tmp38, Scalar(2)) + std::pow(_tmp42, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp44 = _tmp42 * _tmp43;
  const Scalar _tmp45 = _tmp38 * _tmp43;
  const Scalar _tmp46 = _tmp11 + _tmp39;
  const Scalar _tmp47 = _tmp46 + p_init0;
  const Scalar _tmp48 = _tmp47 - p_c(0, 0);
  const Scalar _tmp49 = _tmp24 + _tmp35;
  const Scalar _tmp50 = _tmp49 + p_init1;
  const Scalar _tmp51 = _tmp50 - p_c(1, 0);
  const Scalar _tmp52 =
      std::sqrt(Scalar(std::pow(_tmp48, Scalar(2)) + std::pow(_tmp51, Scalar(2))));
  const Scalar _tmp53 = Scalar(1.0) / (_tmp52);
  const Scalar _tmp54 = Scalar(1.0) / (_tmp48);
  const Scalar _tmp55 = _tmp52 * _tmp54;
  const Scalar _tmp56 = _tmp55 * (_tmp46 * _tmp51 * _tmp53 - _tmp48 * _tmp49 * _tmp53);
  const Scalar _tmp57 = _tmp24 + _tmp32;
  const Scalar _tmp58 = _tmp57 + p_init1;
  const Scalar _tmp59 = _tmp58 - p_d(1, 0);
  const Scalar _tmp60 = _tmp11 + _tmp20;
  const Scalar _tmp61 = _tmp60 + p_init0;
  const Scalar _tmp62 = _tmp61 - p_d(0, 0);
  const Scalar _tmp63 = std::pow(Scalar(std::pow(_tmp59, Scalar(2)) + std::pow(_tmp62, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp64 = _tmp62 * _tmp63;
  const Scalar _tmp65 = _tmp59 * _tmp63;
  const Scalar _tmp66 = _tmp56 * _tmp64 + _tmp57 * _tmp64 - _tmp60 * _tmp65;
  const Scalar _tmp67 = _tmp51 * _tmp54;
  const Scalar _tmp68 = Scalar(1.0) / (_tmp64 * _tmp67 - _tmp65);
  const Scalar _tmp69 = _tmp44 * _tmp67 - _tmp45;
  const Scalar _tmp70 = _tmp68 * _tmp69;
  const Scalar _tmp71 = _tmp36 * _tmp44 - _tmp40 * _tmp45 + _tmp44 * _tmp56 - _tmp66 * _tmp70;
  const Scalar _tmp72 = Scalar(1.0) / (_tmp71);
  const Scalar _tmp73 = Scalar(1.0) * _tmp72;
  const Scalar _tmp74 = Scalar(1.0) * _tmp68;
  const Scalar _tmp75 = _tmp69 * _tmp72 * _tmp74;
  const Scalar _tmp76 = _tmp22 - p_a(0, 0);
  const Scalar _tmp77 = _tmp34 - p_a(1, 0);
  const Scalar _tmp78 = std::pow(Scalar(std::pow(_tmp76, Scalar(2)) + std::pow(_tmp77, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp79 = _tmp77 * _tmp78;
  const Scalar _tmp80 = _tmp76 * _tmp78;
  const Scalar _tmp81 = fh1 * (_tmp21 * _tmp79 - _tmp33 * _tmp80);
  const Scalar _tmp82 = Scalar(0.20999999999999999) * _tmp26 + Scalar(0.20999999999999999) * _tmp27;
  const Scalar _tmp83 = -Scalar(0.010999999999999999) * _tmp17 -
                        Scalar(0.010999999999999999) * _tmp23 + Scalar(-0.010999999999999999);
  const Scalar _tmp84 = Scalar(0.20999999999999999) * _tmp13 - Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp85 = _tmp83 + _tmp84;
  const Scalar _tmp86 = _tmp82 + _tmp85;
  const Scalar _tmp87 = _tmp44 * _tmp86;
  const Scalar _tmp88 = -_tmp82;
  const Scalar _tmp89 = _tmp85 + _tmp88;
  const Scalar _tmp90 = _tmp64 * _tmp86;
  const Scalar _tmp91 = _tmp83 - _tmp84;
  const Scalar _tmp92 = _tmp82 + _tmp91;
  const Scalar _tmp93 = -_tmp64 * _tmp92 + _tmp90;
  const Scalar _tmp94 = _tmp65 * _tmp92 - _tmp67 * _tmp90;
  const Scalar _tmp95 = _tmp45 * _tmp89 - _tmp67 * _tmp87 - _tmp70 * _tmp94;
  const Scalar _tmp96 = Scalar(1.0) * _tmp49;
  const Scalar _tmp97 = -_tmp96;
  const Scalar _tmp98 = Scalar(1.0) / (_tmp57 + _tmp97);
  const Scalar _tmp99 = Scalar(1.0) * _tmp46;
  const Scalar _tmp100 = _tmp98 * (-_tmp60 + _tmp99);
  const Scalar _tmp101 = -_tmp100 * _tmp95 - _tmp44 * _tmp89 - _tmp70 * _tmp93 + _tmp87;
  const Scalar _tmp102 = Scalar(1.0) / (_tmp101);
  const Scalar _tmp103 = _tmp67 * _tmp68;
  const Scalar _tmp104 = _tmp103 * _tmp94 + _tmp67 * _tmp86;
  const Scalar _tmp105 = -_tmp100 * _tmp104 + _tmp103 * _tmp93 - _tmp86;
  const Scalar _tmp106 = _tmp102 * _tmp71;
  const Scalar _tmp107 = _tmp101 * _tmp72;
  const Scalar _tmp108 = _tmp107 * (_tmp103 * _tmp66 - _tmp105 * _tmp106 - _tmp56);
  const Scalar _tmp109 = _tmp102 * (_tmp105 + _tmp108);
  const Scalar _tmp110 = -_tmp109 * _tmp69 - _tmp67;
  const Scalar _tmp111 = _tmp64 * _tmp68;
  const Scalar _tmp112 = _tmp80 * fh1;
  const Scalar _tmp113 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp114 = _tmp100 * _tmp96 + _tmp99;
  const Scalar _tmp115 = 0;
  const Scalar _tmp116 = _tmp115 * _tmp70;
  const Scalar _tmp117 = _tmp55 * (_tmp115 * _tmp44 - _tmp116 * _tmp64);
  const Scalar _tmp118 = _tmp74 * _tmp94;
  const Scalar _tmp119 = _tmp100 * _tmp118 - _tmp74 * _tmp93;
  const Scalar _tmp120 = _tmp107 * (-_tmp106 * _tmp119 - _tmp66 * _tmp74);
  const Scalar _tmp121 = _tmp102 * (_tmp119 + _tmp120);
  const Scalar _tmp122 = -_tmp121 * _tmp69 + Scalar(1.0);
  const Scalar _tmp123 = _tmp79 * fh1;
  const Scalar _tmp124 = -_tmp112 * _tmp55 * (_tmp109 * _tmp44 + _tmp110 * _tmp111 + Scalar(1.0)) -
                         _tmp113 * _tmp117 -
                         _tmp123 * _tmp55 * (_tmp111 * _tmp122 + _tmp121 * _tmp44) -
                         _tmp55 * _tmp81 * (_tmp44 * _tmp73 - _tmp64 * _tmp75);
  const Scalar _tmp125 = Scalar(1.0) / (_tmp124);
  const Scalar _tmp126 = _tmp36 + _tmp97;
  const Scalar _tmp127 = _tmp100 * _tmp126;
  const Scalar _tmp128 = Scalar(1.0) / (-_tmp127 - _tmp40 + _tmp99);
  const Scalar _tmp129 = _tmp114 * _tmp128;
  const Scalar _tmp130 = _tmp98 * (-_tmp115 * _tmp95 - _tmp126 * _tmp129 + _tmp97);
  const Scalar _tmp131 = -Scalar(1.0) * _tmp129 - Scalar(1.0) * _tmp130 + Scalar(1.0);
  const Scalar _tmp132 = fh1 * (_tmp88 + _tmp91);
  const Scalar _tmp133 = _tmp132 * _tmp80 + Scalar(5.1796800000000003) * _tmp15 + _tmp21 * fv1;
  const Scalar _tmp134 = Scalar(1.0) * _tmp128;
  const Scalar _tmp135 = _tmp126 * _tmp98;
  const Scalar _tmp136 = Scalar(1.0) * _tmp134 * _tmp135 - Scalar(1.0) * _tmp134;
  const Scalar _tmp137 = -_tmp132 * _tmp79 - Scalar(5.1796800000000003) * _tmp28 - _tmp33 * fv1;
  const Scalar _tmp138 = _tmp127 * _tmp134 + Scalar(1.0);
  const Scalar _tmp139 = Scalar(1.0) * _tmp98;
  const Scalar _tmp140 = _tmp100 * _tmp134;
  const Scalar _tmp141 = -Scalar(1.0) * _tmp138 * _tmp139 + Scalar(1.0) * _tmp140;
  const Scalar _tmp142 = _tmp126 * _tmp128;
  const Scalar _tmp143 = -_tmp118 + _tmp120 * _tmp142 - _tmp121 * _tmp95;
  const Scalar _tmp144 = _tmp104 + _tmp108 * _tmp142 - _tmp109 * _tmp95;
  const Scalar _tmp145 = _tmp107 * _tmp134;
  const Scalar _tmp146 = _tmp126 * _tmp145 - _tmp73 * _tmp95;
  const Scalar _tmp147 =
      Scalar(1.0) * _tmp112 * (_tmp108 * _tmp134 - _tmp139 * _tmp144) + _tmp113 * _tmp131 +
      Scalar(1.0) * _tmp123 * (_tmp120 * _tmp134 - _tmp139 * _tmp143) + _tmp133 * _tmp136 +
      _tmp137 * _tmp141 + Scalar(1.0) * _tmp81 * (-_tmp139 * _tmp146 + _tmp145);
  const Scalar _tmp148 = std::asinh(_tmp125 * _tmp147);
  const Scalar _tmp149 = Scalar(1.0) * _tmp148;
  const Scalar _tmp150 = Scalar(9.6622558468725703) * _tmp124;
  const Scalar _tmp151 =
      -_tmp148 * _tmp150 - std::sqrt(Scalar(std::pow(Scalar(-_tmp47 + p_c(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp50 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp152 = Scalar(0.1034955) * _tmp125;
  const Scalar _tmp153 = _tmp151 * _tmp152;
  const Scalar _tmp154 = Scalar(9.6622558468725703) * _tmp117;
  const Scalar _tmp155 = std::pow(_tmp124, Scalar(-2));
  const Scalar _tmp156 = _tmp117 * _tmp155;
  const Scalar _tmp157 = Scalar(0.1034955) * _tmp156;
  const Scalar _tmp158 = _tmp24 + _tmp29 + _tmp31;
  const Scalar _tmp159 =
      (_tmp125 * (-_tmp131 + _tmp136 * _tmp21 + _tmp141 * _tmp158) - _tmp147 * _tmp156) /
      std::sqrt(Scalar(std::pow(_tmp147, Scalar(2)) * _tmp155 + 1));
  const Scalar _tmp160 = _tmp113 * _tmp115;
  const Scalar _tmp161 =
      _tmp110 * _tmp112 * _tmp68 + _tmp122 * _tmp123 * _tmp68 - _tmp160 * _tmp70 - _tmp75 * _tmp81;
  const Scalar _tmp162 = Scalar(1.0) / (_tmp161);
  const Scalar _tmp163 = _tmp133 * _tmp134;
  const Scalar _tmp164 = _tmp138 * _tmp98;
  const Scalar _tmp165 = _tmp112 * _tmp144 * _tmp98 + _tmp113 * _tmp130 +
                         _tmp123 * _tmp143 * _tmp98 - _tmp135 * _tmp163 + _tmp137 * _tmp164 +
                         _tmp146 * _tmp81 * _tmp98;
  const Scalar _tmp166 = std::asinh(_tmp162 * _tmp165);
  const Scalar _tmp167 = Scalar(1.0) * _tmp166;
  const Scalar _tmp168 = std::pow(_tmp161, Scalar(-2));
  const Scalar _tmp169 = _tmp116 * _tmp168;
  const Scalar _tmp170 = _tmp134 * _tmp21;
  const Scalar _tmp171 =
      (_tmp162 * (-_tmp130 - _tmp135 * _tmp170 + _tmp158 * _tmp164) - _tmp165 * _tmp169) /
      std::sqrt(Scalar(std::pow(_tmp165, Scalar(2)) * _tmp168 + 1));
  const Scalar _tmp172 = Scalar(9.6622558468725703) * _tmp161;
  const Scalar _tmp173 =
      -_tmp166 * _tmp172 - std::sqrt(Scalar(std::pow(Scalar(-_tmp58 + p_d(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp61 + p_d(0, 0)), Scalar(2))));
  const Scalar _tmp174 = Scalar(0.1034955) * _tmp169;
  const Scalar _tmp175 = Scalar(9.6622558468725703) * _tmp115;
  const Scalar _tmp176 = _tmp175 * _tmp70;
  const Scalar _tmp177 = Scalar(0.1034955) * _tmp162;
  const Scalar _tmp178 = _tmp173 * _tmp177;
  const Scalar _tmp179 = _tmp109 * _tmp112 + _tmp121 * _tmp123 + _tmp160 + _tmp73 * _tmp81;
  const Scalar _tmp180 = Scalar(1.0) / (_tmp179);
  const Scalar _tmp181 = -_tmp108 * _tmp112 * _tmp128 + _tmp113 * _tmp129 -
                         _tmp120 * _tmp123 * _tmp128 - _tmp137 * _tmp140 - _tmp145 * _tmp81 +
                         _tmp163;
  const Scalar _tmp182 = std::asinh(_tmp180 * _tmp181);
  const Scalar _tmp183 = Scalar(1.0) * _tmp182;
  const Scalar _tmp184 = std::pow(_tmp179, Scalar(-2));
  const Scalar _tmp185 = _tmp115 * _tmp184;
  const Scalar _tmp186 = (_tmp180 * (-_tmp129 - _tmp140 * _tmp158 + _tmp170) + _tmp181 * _tmp185) /
                         std::sqrt(Scalar(std::pow(_tmp181, Scalar(2)) * _tmp184 + 1));
  const Scalar _tmp187 = Scalar(9.6622558468725703) * _tmp179;
  const Scalar _tmp188 = Scalar(0.1034955) * _tmp180;
  const Scalar _tmp189 =
      -_tmp182 * _tmp187 - std::sqrt(Scalar(std::pow(Scalar(-_tmp37 + p_b(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp41 + p_b(0, 0)), Scalar(2))));
  const Scalar _tmp190 = Scalar(0.1034955) * _tmp185;
  const Scalar _tmp191 = _tmp188 * _tmp189;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      -_tmp3 *
      (_tmp2 * std::sinh(Scalar(1.0) * _tmp1) +
       _tmp2 * std::sinh(Scalar(0.1034955) * _tmp0 *
                         (-_tmp1 * _tmp3 -
                          std::sqrt(Scalar(std::pow(Scalar(-_tmp22 + p_a(0, 0)), Scalar(2)) +
                                           std::pow(Scalar(-_tmp34 + p_a(1, 0)), Scalar(2)))))));
  _res(1, 0) =
      -_tmp150 * (-_tmp157 * p_c(2, 0) + Scalar(1.0) * _tmp159 * std::sinh(_tmp149) -
                  (-_tmp151 * _tmp157 + _tmp152 * (-_tmp148 * _tmp154 - _tmp150 * _tmp159)) *
                      std::sinh(_tmp153)) -
      _tmp154 * (_tmp152 * p_c(2, 0) + std::cosh(_tmp149) - std::cosh(_tmp153));
  _res(2, 0) =
      -_tmp172 * (Scalar(1.0) * _tmp171 * std::sinh(_tmp167) - _tmp174 * p_d(2, 0) -
                  (-_tmp173 * _tmp174 + _tmp177 * (-_tmp166 * _tmp176 - _tmp171 * _tmp172)) *
                      std::sinh(_tmp178)) -
      _tmp176 * (_tmp177 * p_d(2, 0) + std::cosh(_tmp167) - std::cosh(_tmp178));
  _res(3, 0) = _tmp175 * (_tmp188 * p_b(2, 0) + std::cosh(_tmp183) - std::cosh(_tmp191)) -
               _tmp187 * (Scalar(1.0) * _tmp186 * std::sinh(_tmp183) + _tmp190 * p_b(2, 0) -
                          (_tmp188 * (_tmp175 * _tmp182 - _tmp186 * _tmp187) + _tmp189 * _tmp190) *
                              std::sinh(_tmp191));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
