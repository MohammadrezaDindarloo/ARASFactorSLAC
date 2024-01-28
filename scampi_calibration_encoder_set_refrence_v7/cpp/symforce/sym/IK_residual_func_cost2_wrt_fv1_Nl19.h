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
 * Symbolic function: IK_residual_func_cost2_wrt_fv1_Nl19
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2WrtFv1Nl19(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const Scalar p_init0, const Scalar p_init1, const Scalar p_init2, const Scalar rot_init_x,
    const Scalar rot_init_y, const Scalar rot_init_z, const Scalar rot_init_w,
    const Scalar epsilon) {
  // Total ops: 578

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
  const Scalar _tmp3 = _DeltaRot[0] * rot_init_z + _DeltaRot[1] * rot_init_w -
                       _DeltaRot[2] * rot_init_x + _DeltaRot[3] * rot_init_y;
  const Scalar _tmp4 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = -_DeltaRot[0] * rot_init_y + _DeltaRot[1] * rot_init_x +
                       _DeltaRot[2] * rot_init_w + _DeltaRot[3] * rot_init_z;
  const Scalar _tmp6 = 1 - 2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp7 = Scalar(0.20999999999999999) * _tmp4 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp8 = -_tmp7;
  const Scalar _tmp9 = _DeltaRot[0] * rot_init_w - _DeltaRot[1] * rot_init_z +
                       _DeltaRot[2] * rot_init_y + _DeltaRot[3] * rot_init_x;
  const Scalar _tmp10 = 2 * _tmp9;
  const Scalar _tmp11 = _tmp10 * _tmp5;
  const Scalar _tmp12 = -2 * _DeltaRot[0] * rot_init_x - 2 * _DeltaRot[1] * rot_init_y -
                        2 * _DeltaRot[2] * rot_init_z + 2 * _DeltaRot[3] * rot_init_w;
  const Scalar _tmp13 = _tmp12 * _tmp3;
  const Scalar _tmp14 = _tmp11 + _tmp13;
  const Scalar _tmp15 = -Scalar(0.010999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp10 * _tmp3;
  const Scalar _tmp17 = _tmp12 * _tmp5;
  const Scalar _tmp18 = Scalar(0.20999999999999999) * _tmp16 - Scalar(0.20999999999999999) * _tmp17;
  const Scalar _tmp19 = _tmp15 + _tmp18;
  const Scalar _tmp20 = _tmp19 + _tmp8;
  const Scalar _tmp21 = _tmp20 + p_init0;
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp16 + Scalar(0.20999999999999999) * _tmp17;
  const Scalar _tmp23 = -_tmp22;
  const Scalar _tmp24 = -2 * std::pow(_tmp9, Scalar(2));
  const Scalar _tmp25 = Scalar(0.20999999999999999) * _tmp24 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp26 = 2 * _tmp3 * _tmp5;
  const Scalar _tmp27 = _tmp12 * _tmp9;
  const Scalar _tmp28 = _tmp26 - _tmp27;
  const Scalar _tmp29 = Scalar(0.010999999999999999) * _tmp28;
  const Scalar _tmp30 = -_tmp29;
  const Scalar _tmp31 = _tmp25 + _tmp30;
  const Scalar _tmp32 = _tmp23 + _tmp31;
  const Scalar _tmp33 = _tmp32 + p_init1;
  const Scalar _tmp34 = Scalar(9.6622558468725703) * fh1;
  const Scalar _tmp35 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp36 = Scalar(0.20999999999999999) * _tmp11 - Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp37 = -_tmp36;
  const Scalar _tmp38 = -Scalar(0.010999999999999999) * _tmp24 -
                        Scalar(0.010999999999999999) * _tmp4 + Scalar(-0.010999999999999999);
  const Scalar _tmp39 = Scalar(0.20999999999999999) * _tmp26 + Scalar(0.20999999999999999) * _tmp27;
  const Scalar _tmp40 = _tmp38 - _tmp39;
  const Scalar _tmp41 = _tmp37 + _tmp40;
  const Scalar _tmp42 = -_tmp25;
  const Scalar _tmp43 = _tmp22 + _tmp42;
  const Scalar _tmp44 = _tmp30 + _tmp43;
  const Scalar _tmp45 = _tmp44 + p_init1;
  const Scalar _tmp46 = _tmp45 - p_b(1, 0);
  const Scalar _tmp47 = _tmp15 - _tmp18;
  const Scalar _tmp48 = _tmp47 + _tmp7;
  const Scalar _tmp49 = _tmp48 + p_init0;
  const Scalar _tmp50 = _tmp49 - p_b(0, 0);
  const Scalar _tmp51 = std::pow(Scalar(std::pow(_tmp46, Scalar(2)) + std::pow(_tmp50, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp52 = _tmp50 * _tmp51;
  const Scalar _tmp53 = _tmp36 + _tmp40;
  const Scalar _tmp54 = _tmp23 + _tmp30 + _tmp42;
  const Scalar _tmp55 = _tmp54 + p_init1;
  const Scalar _tmp56 = _tmp55 - p_a(1, 0);
  const Scalar _tmp57 = _tmp47 + _tmp8;
  const Scalar _tmp58 = _tmp57 + p_init0;
  const Scalar _tmp59 = _tmp58 - p_a(0, 0);
  const Scalar _tmp60 = Scalar(1.0) / (_tmp59);
  const Scalar _tmp61 = _tmp56 * _tmp60;
  const Scalar _tmp62 = _tmp41 * _tmp61;
  const Scalar _tmp63 = _tmp46 * _tmp51;
  const Scalar _tmp64 = _tmp19 + _tmp7;
  const Scalar _tmp65 = _tmp64 + p_init0;
  const Scalar _tmp66 = _tmp65 - p_c(0, 0);
  const Scalar _tmp67 = _tmp22 + _tmp31;
  const Scalar _tmp68 = _tmp67 + p_init1;
  const Scalar _tmp69 = _tmp68 - p_c(1, 0);
  const Scalar _tmp70 = std::pow(Scalar(std::pow(_tmp66, Scalar(2)) + std::pow(_tmp69, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp71 = _tmp66 * _tmp70;
  const Scalar _tmp72 = _tmp41 * _tmp71;
  const Scalar _tmp73 = _tmp38 + _tmp39;
  const Scalar _tmp74 = _tmp36 + _tmp73;
  const Scalar _tmp75 = _tmp69 * _tmp70;
  const Scalar _tmp76 = -_tmp61 * _tmp72 + _tmp74 * _tmp75;
  const Scalar _tmp77 = Scalar(1.0) / (_tmp61 * _tmp71 - _tmp75);
  const Scalar _tmp78 = _tmp52 * _tmp61 - _tmp63;
  const Scalar _tmp79 = _tmp77 * _tmp78;
  const Scalar _tmp80 = -_tmp52 * _tmp62 + _tmp53 * _tmp63 - _tmp76 * _tmp79;
  const Scalar _tmp81 = Scalar(1.0) * _tmp54;
  const Scalar _tmp82 = -_tmp81;
  const Scalar _tmp83 = Scalar(1.0) / (_tmp67 + _tmp82);
  const Scalar _tmp84 = Scalar(1.0) * _tmp57;
  const Scalar _tmp85 = -_tmp64 + _tmp84;
  const Scalar _tmp86 = _tmp83 * _tmp85;
  const Scalar _tmp87 = -_tmp71 * _tmp74 + _tmp72;
  const Scalar _tmp88 = _tmp41 * _tmp52 - _tmp52 * _tmp53 - _tmp79 * _tmp87 - _tmp80 * _tmp86;
  const Scalar _tmp89 = Scalar(1.0) / (_tmp88);
  const Scalar _tmp90 = _tmp81 * _tmp86 + _tmp84;
  const Scalar _tmp91 = 0;
  const Scalar _tmp92 = _tmp79 * _tmp91;
  const Scalar _tmp93 =
      std::sqrt(Scalar(std::pow(_tmp56, Scalar(2)) + std::pow(_tmp59, Scalar(2))));
  const Scalar _tmp94 = _tmp60 * _tmp93;
  const Scalar _tmp95 = _tmp94 * (_tmp52 * _tmp91 - _tmp71 * _tmp92);
  const Scalar _tmp96 = Scalar(1.0) / (_tmp93);
  const Scalar _tmp97 = _tmp94 * (-_tmp54 * _tmp59 * _tmp96 + _tmp56 * _tmp57 * _tmp96);
  const Scalar _tmp98 = -_tmp64 * _tmp75 + _tmp67 * _tmp71 + _tmp71 * _tmp97;
  const Scalar _tmp99 = _tmp44 * _tmp52 - _tmp48 * _tmp63 + _tmp52 * _tmp97 - _tmp79 * _tmp98;
  const Scalar _tmp100 = Scalar(1.0) / (_tmp99);
  const Scalar _tmp101 = Scalar(1.0) * _tmp100;
  const Scalar _tmp102 = _tmp33 - p_d(1, 0);
  const Scalar _tmp103 = _tmp21 - p_d(0, 0);
  const Scalar _tmp104 =
      std::pow(Scalar(std::pow(_tmp102, Scalar(2)) + std::pow(_tmp103, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp105 = _tmp103 * _tmp104;
  const Scalar _tmp106 = _tmp102 * _tmp104;
  const Scalar _tmp107 = fh1 * (-_tmp105 * _tmp32 + _tmp106 * _tmp20);
  const Scalar _tmp108 = Scalar(1.0) * _tmp83;
  const Scalar _tmp109 = Scalar(1.0) * _tmp77;
  const Scalar _tmp110 = _tmp108 * _tmp76 * _tmp77 * _tmp85 - _tmp109 * _tmp87;
  const Scalar _tmp111 = _tmp89 * _tmp99;
  const Scalar _tmp112 = _tmp100 * _tmp88;
  const Scalar _tmp113 = _tmp112 * (-_tmp109 * _tmp98 - _tmp110 * _tmp111);
  const Scalar _tmp114 = _tmp89 * (_tmp110 + _tmp113);
  const Scalar _tmp115 = -_tmp114 * _tmp78 + Scalar(1.0);
  const Scalar _tmp116 = _tmp71 * _tmp77;
  const Scalar _tmp117 = _tmp106 * fh1;
  const Scalar _tmp118 = _tmp61 * _tmp77;
  const Scalar _tmp119 = _tmp118 * _tmp76 + _tmp62;
  const Scalar _tmp120 = _tmp118 * _tmp87 - _tmp119 * _tmp86 - _tmp41;
  const Scalar _tmp121 = _tmp112 * (-_tmp111 * _tmp120 + _tmp118 * _tmp98 - _tmp97);
  const Scalar _tmp122 = _tmp89 * (_tmp120 + _tmp121);
  const Scalar _tmp123 = -_tmp122 * _tmp78 - _tmp61;
  const Scalar _tmp124 = _tmp105 * fh1;
  const Scalar _tmp125 = -_tmp107 * _tmp94 * (_tmp101 * _tmp52 - _tmp101 * _tmp71 * _tmp79) -
                         _tmp117 * _tmp94 * (_tmp114 * _tmp52 + _tmp115 * _tmp116) -
                         _tmp124 * _tmp94 * (_tmp116 * _tmp123 + _tmp122 * _tmp52 + Scalar(1.0)) -
                         _tmp35 * _tmp95;
  const Scalar _tmp126 = Scalar(1.0) / (_tmp125);
  const Scalar _tmp127 = _tmp44 + _tmp82;
  const Scalar _tmp128 = _tmp127 * _tmp86;
  const Scalar _tmp129 = Scalar(1.0) / (-_tmp128 - _tmp48 + _tmp84);
  const Scalar _tmp130 = Scalar(1.0) * _tmp129;
  const Scalar _tmp131 = _tmp127 * _tmp129;
  const Scalar _tmp132 = -_tmp109 * _tmp76 + _tmp113 * _tmp131 - _tmp114 * _tmp80;
  const Scalar _tmp133 = Scalar(1.0) * fh1;
  const Scalar _tmp134 = _tmp129 * _tmp90;
  const Scalar _tmp135 = _tmp83 * (-_tmp127 * _tmp134 - _tmp80 * _tmp91 + _tmp82);
  const Scalar _tmp136 = -Scalar(1.0) * _tmp134 - Scalar(1.0) * _tmp135 + Scalar(1.0);
  const Scalar _tmp137 = fh1 * (_tmp37 + _tmp73);
  const Scalar _tmp138 = _tmp105 * _tmp137 + Scalar(5.1796800000000003) * _tmp14 + _tmp20 * fv1;
  const Scalar _tmp139 = _tmp127 * _tmp83;
  const Scalar _tmp140 = Scalar(1.0) * _tmp130 * _tmp139 - Scalar(1.0) * _tmp130;
  const Scalar _tmp141 = _tmp112 * _tmp130;
  const Scalar _tmp142 = -_tmp101 * _tmp80 + _tmp127 * _tmp141;
  const Scalar _tmp143 = _tmp119 + _tmp121 * _tmp131 - _tmp122 * _tmp80;
  const Scalar _tmp144 = -_tmp106 * _tmp137 - Scalar(5.1796800000000003) * _tmp28 - _tmp32 * fv1;
  const Scalar _tmp145 = _tmp130 * _tmp86;
  const Scalar _tmp146 = _tmp128 * _tmp130 + Scalar(1.0);
  const Scalar _tmp147 = -Scalar(1.0) * _tmp108 * _tmp146 + Scalar(1.0) * _tmp145;
  const Scalar _tmp148 = _tmp105 * _tmp133 * (-_tmp108 * _tmp143 + _tmp121 * _tmp130) +
                         _tmp106 * _tmp133 * (-_tmp108 * _tmp132 + _tmp113 * _tmp130) +
                         Scalar(1.0) * _tmp107 * (-_tmp108 * _tmp142 + _tmp141) + _tmp136 * _tmp35 +
                         _tmp138 * _tmp140 + _tmp144 * _tmp147;
  const Scalar _tmp149 = std::asinh(_tmp126 * _tmp148);
  const Scalar _tmp150 = Scalar(9.6622558468725703) * _tmp125;
  const Scalar _tmp151 =
      -_tmp149 * _tmp150 - std::sqrt(Scalar(std::pow(Scalar(-_tmp55 + p_a(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp58 + p_a(0, 0)), Scalar(2))));
  const Scalar _tmp152 = std::pow(_tmp125, Scalar(-2));
  const Scalar _tmp153 = _tmp152 * _tmp95;
  const Scalar _tmp154 = Scalar(9.6622558468725703) * _tmp95;
  const Scalar _tmp155 = _tmp29 + _tmp43;
  const Scalar _tmp156 =
      (_tmp126 * (-_tmp136 + _tmp140 * _tmp20 + _tmp147 * _tmp155) - _tmp148 * _tmp153) /
      std::sqrt(Scalar(std::pow(_tmp148, Scalar(2)) * _tmp152 + 1));
  const Scalar _tmp157 = Scalar(0.1034955) * _tmp126;
  const Scalar _tmp158 = _tmp151 * _tmp157;
  const Scalar _tmp159 = Scalar(1.0) * _tmp149;
  const Scalar _tmp160 = _tmp101 * _tmp107;
  const Scalar _tmp161 = _tmp35 * _tmp91;
  const Scalar _tmp162 =
      _tmp115 * _tmp117 * _tmp77 + _tmp123 * _tmp124 * _tmp77 - _tmp160 * _tmp79 - _tmp161 * _tmp79;
  const Scalar _tmp163 = Scalar(1.0) / (_tmp162);
  const Scalar _tmp164 = _tmp146 * _tmp83;
  const Scalar _tmp165 = _tmp83 * fh1;
  const Scalar _tmp166 = _tmp130 * _tmp138;
  const Scalar _tmp167 = _tmp105 * _tmp143 * _tmp165 + _tmp106 * _tmp132 * _tmp165 +
                         _tmp107 * _tmp142 * _tmp83 + _tmp135 * _tmp35 - _tmp139 * _tmp166 +
                         _tmp144 * _tmp164;
  const Scalar _tmp168 = std::asinh(_tmp163 * _tmp167);
  const Scalar _tmp169 = Scalar(9.6622558468725703) * _tmp162;
  const Scalar _tmp170 =
      -_tmp168 * _tmp169 - std::sqrt(Scalar(std::pow(Scalar(-_tmp65 + p_c(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp68 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp171 = std::pow(_tmp162, Scalar(-2));
  const Scalar _tmp172 = _tmp171 * _tmp92;
  const Scalar _tmp173 = Scalar(9.6622558468725703) * _tmp91;
  const Scalar _tmp174 = _tmp173 * _tmp79;
  const Scalar _tmp175 = _tmp130 * _tmp20;
  const Scalar _tmp176 =
      (_tmp163 * (-_tmp135 - _tmp139 * _tmp175 + _tmp155 * _tmp164) - _tmp167 * _tmp172) /
      std::sqrt(Scalar(std::pow(_tmp167, Scalar(2)) * _tmp171 + 1));
  const Scalar _tmp177 = Scalar(0.1034955) * _tmp163;
  const Scalar _tmp178 = _tmp170 * _tmp177;
  const Scalar _tmp179 = Scalar(1.0) * _tmp168;
  const Scalar _tmp180 = -_tmp107 * _tmp141 - _tmp113 * _tmp117 * _tmp129 -
                         _tmp121 * _tmp124 * _tmp129 + _tmp134 * _tmp35 - _tmp144 * _tmp145 +
                         _tmp166;
  const Scalar _tmp181 = _tmp114 * _tmp117 + _tmp122 * _tmp124 + _tmp160 + _tmp161;
  const Scalar _tmp182 = Scalar(1.0) / (_tmp181);
  const Scalar _tmp183 = std::asinh(_tmp180 * _tmp182);
  const Scalar _tmp184 = Scalar(1.0) * _tmp183;
  const Scalar _tmp185 = std::pow(_tmp181, Scalar(-2));
  const Scalar _tmp186 = _tmp185 * _tmp91;
  const Scalar _tmp187 = (_tmp180 * _tmp186 + _tmp182 * (-_tmp134 - _tmp145 * _tmp155 + _tmp175)) /
                         std::sqrt(Scalar(std::pow(_tmp180, Scalar(2)) * _tmp185 + 1));
  const Scalar _tmp188 = Scalar(9.6622558468725703) * _tmp181;
  const Scalar _tmp189 =
      -_tmp183 * _tmp188 - std::sqrt(Scalar(std::pow(Scalar(-_tmp45 + p_b(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp49 + p_b(0, 0)), Scalar(2))));
  const Scalar _tmp190 = Scalar(0.1034955) * _tmp182;
  const Scalar _tmp191 = _tmp189 * _tmp190;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      _tmp34 *
      (-_tmp2 * std::cosh(Scalar(1.0) * _tmp1) +
       _tmp2 * std::cosh(Scalar(0.1034955) * _tmp0 *
                         (-_tmp1 * _tmp34 -
                          std::sqrt(Scalar(std::pow(Scalar(-_tmp21 + p_d(0, 0)), Scalar(2)) +
                                           std::pow(Scalar(-_tmp33 + p_d(1, 0)), Scalar(2)))))));
  _res(1, 0) = _tmp150 * (-Scalar(1.0) * _tmp156 * std::cosh(_tmp159) -
                          (-Scalar(0.1034955) * _tmp151 * _tmp153 +
                           _tmp157 * (-_tmp149 * _tmp154 - _tmp150 * _tmp156)) *
                              std::cosh(_tmp158)) +
               _tmp154 * (-std::sinh(_tmp158) - std::sinh(_tmp159));
  _res(2, 0) = _tmp169 * (-Scalar(1.0) * _tmp176 * std::cosh(_tmp179) -
                          (-Scalar(0.1034955) * _tmp170 * _tmp172 +
                           _tmp177 * (-_tmp168 * _tmp174 - _tmp169 * _tmp176)) *
                              std::cosh(_tmp178)) +
               _tmp174 * (-std::sinh(_tmp178) - std::sinh(_tmp179));
  _res(3, 0) = -_tmp173 * (-std::sinh(_tmp184) - std::sinh(_tmp191)) +
               _tmp188 * (-Scalar(1.0) * _tmp187 * std::cosh(_tmp184) -
                          (Scalar(0.1034955) * _tmp186 * _tmp189 +
                           _tmp190 * (_tmp173 * _tmp183 - _tmp187 * _tmp188)) *
                              std::cosh(_tmp191));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
