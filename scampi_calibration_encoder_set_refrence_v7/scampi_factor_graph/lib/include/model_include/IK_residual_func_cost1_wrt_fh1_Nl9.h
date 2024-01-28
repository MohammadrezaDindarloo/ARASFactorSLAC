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
 * Symbolic function: IK_residual_func_cost1_wrt_fh1_Nl9
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1WrtFh1Nl9(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const Scalar p_init0, const Scalar p_init1, const Scalar p_init2, const Scalar rot_init_x,
    const Scalar rot_init_y, const Scalar rot_init_z, const Scalar rot_init_w,
    const Scalar epsilon) {
  // Total ops: 640

  // Unused inputs
  (void)p_init2;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();

  // Intermediate terms (217)
  const Scalar _tmp0 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp1 = Scalar(0.1034955) * p_b(2, 0);
  const Scalar _tmp2 = _DeltaRot[0] * rot_init_w - _DeltaRot[1] * rot_init_z +
                       _DeltaRot[2] * rot_init_y + _DeltaRot[3] * rot_init_x;
  const Scalar _tmp3 = _DeltaRot[0] * rot_init_z + _DeltaRot[1] * rot_init_w -
                       _DeltaRot[2] * rot_init_x + _DeltaRot[3] * rot_init_y;
  const Scalar _tmp4 = 2 * _tmp3;
  const Scalar _tmp5 = _tmp2 * _tmp4;
  const Scalar _tmp6 = -_DeltaRot[0] * rot_init_x - _DeltaRot[1] * rot_init_y -
                       _DeltaRot[2] * rot_init_z + _DeltaRot[3] * rot_init_w;
  const Scalar _tmp7 = -_DeltaRot[0] * rot_init_y + _DeltaRot[1] * rot_init_x +
                       _DeltaRot[2] * rot_init_w + _DeltaRot[3] * rot_init_z;
  const Scalar _tmp8 = 2 * _tmp7;
  const Scalar _tmp9 = _tmp6 * _tmp8;
  const Scalar _tmp10 = Scalar(0.20999999999999999) * _tmp5 - Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp11 = -_tmp10;
  const Scalar _tmp12 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp13 = 1 - 2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp2 * _tmp8;
  const Scalar _tmp16 = _tmp4 * _tmp6;
  const Scalar _tmp17 = _tmp15 + _tmp16;
  const Scalar _tmp18 = -Scalar(0.010999999999999999) * _tmp17;
  const Scalar _tmp19 = _tmp14 + _tmp18;
  const Scalar _tmp20 = _tmp11 + _tmp19;
  const Scalar _tmp21 = _tmp20 + p_init0;
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp5 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp23 = _tmp3 * _tmp8;
  const Scalar _tmp24 = 2 * _tmp2 * _tmp6;
  const Scalar _tmp25 = _tmp23 - _tmp24;
  const Scalar _tmp26 = -Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = -2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp28 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp27;
  const Scalar _tmp29 = _tmp26 - _tmp28;
  const Scalar _tmp30 = _tmp22 + _tmp29;
  const Scalar _tmp31 = _tmp30 + p_init1;
  const Scalar _tmp32 = Scalar(1.0) / (fh1);
  const Scalar _tmp33 = _tmp32 * fv1;
  const Scalar _tmp34 = std::asinh(_tmp33);
  const Scalar _tmp35 = Scalar(9.6622558468725703) * _tmp34;
  const Scalar _tmp36 =
      -_tmp35 * fh1 - std::sqrt(Scalar(std::pow(Scalar(-_tmp21 + p_b(0, 0)), Scalar(2)) +
                                       std::pow(Scalar(-_tmp31 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp37 = Scalar(0.1034955) * _tmp32;
  const Scalar _tmp38 = _tmp36 * _tmp37;
  const Scalar _tmp39 =
      std::pow(Scalar(_tmp0 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp40 = Scalar(1.0) * _tmp34;
  const Scalar _tmp41 = _tmp26 + _tmp28;
  const Scalar _tmp42 = _tmp22 + _tmp41;
  const Scalar _tmp43 = _tmp42 + p_init1;
  const Scalar _tmp44 = _tmp43 - p_c(1, 0);
  const Scalar _tmp45 = _tmp10 + _tmp19;
  const Scalar _tmp46 = _tmp45 + p_init0;
  const Scalar _tmp47 = _tmp46 - p_c(0, 0);
  const Scalar _tmp48 = Scalar(1.0) / (_tmp47);
  const Scalar _tmp49 = _tmp44 * _tmp48;
  const Scalar _tmp50 = -_tmp22;
  const Scalar _tmp51 = _tmp41 + _tmp50;
  const Scalar _tmp52 = _tmp51 + p_init1;
  const Scalar _tmp53 = _tmp52 - p_d(1, 0);
  const Scalar _tmp54 = -_tmp14 + _tmp18;
  const Scalar _tmp55 = _tmp10 + _tmp54;
  const Scalar _tmp56 = _tmp55 + p_init0;
  const Scalar _tmp57 = _tmp56 - p_d(0, 0);
  const Scalar _tmp58 = std::pow(Scalar(std::pow(_tmp53, Scalar(2)) + std::pow(_tmp57, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp59 = _tmp53 * _tmp58;
  const Scalar _tmp60 = _tmp57 * _tmp58;
  const Scalar _tmp61 = Scalar(1.0) / (_tmp49 * _tmp60 - _tmp59);
  const Scalar _tmp62 =
      std::sqrt(Scalar(std::pow(_tmp44, Scalar(2)) + std::pow(_tmp47, Scalar(2))));
  const Scalar _tmp63 = Scalar(1.0) / (_tmp62);
  const Scalar _tmp64 = _tmp48 * _tmp62;
  const Scalar _tmp65 = _tmp64 * (-_tmp42 * _tmp47 * _tmp63 + _tmp44 * _tmp45 * _tmp63);
  const Scalar _tmp66 = _tmp61 * (_tmp51 * _tmp60 - _tmp55 * _tmp59 + _tmp60 * _tmp65);
  const Scalar _tmp67 = Scalar(0.20999999999999999) * _tmp15 - Scalar(0.20999999999999999) * _tmp16;
  const Scalar _tmp68 = -Scalar(0.010999999999999999) * _tmp12 -
                        Scalar(0.010999999999999999) * _tmp27 + Scalar(-0.010999999999999999);
  const Scalar _tmp69 = Scalar(0.20999999999999999) * _tmp23 + Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp70 = _tmp68 + _tmp69;
  const Scalar _tmp71 = _tmp67 + _tmp70;
  const Scalar _tmp72 = -_tmp67;
  const Scalar _tmp73 = _tmp70 + _tmp72;
  const Scalar _tmp74 = _tmp60 * _tmp71;
  const Scalar _tmp75 = _tmp61 * (-_tmp49 * _tmp74 + _tmp59 * _tmp73);
  const Scalar _tmp76 = _tmp49 * _tmp71 + _tmp49 * _tmp75;
  const Scalar _tmp77 = Scalar(1.0) * _tmp42;
  const Scalar _tmp78 = -_tmp77;
  const Scalar _tmp79 = Scalar(1.0) / (_tmp51 + _tmp78);
  const Scalar _tmp80 = Scalar(1.0) * _tmp45;
  const Scalar _tmp81 = _tmp79 * (-_tmp55 + _tmp80);
  const Scalar _tmp82 = _tmp61 * (-_tmp60 * _tmp73 + _tmp74);
  const Scalar _tmp83 = _tmp49 * _tmp82 - _tmp71 - _tmp76 * _tmp81;
  const Scalar _tmp84 = _tmp11 + _tmp54;
  const Scalar _tmp85 = _tmp84 + p_init0;
  const Scalar _tmp86 = _tmp85 - p_a(0, 0);
  const Scalar _tmp87 = _tmp29 + _tmp50;
  const Scalar _tmp88 = _tmp87 + p_init1;
  const Scalar _tmp89 = _tmp88 - p_a(1, 0);
  const Scalar _tmp90 = std::pow(Scalar(std::pow(_tmp86, Scalar(2)) + std::pow(_tmp89, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp91 = _tmp86 * _tmp90;
  const Scalar _tmp92 = _tmp89 * _tmp90;
  const Scalar _tmp93 = _tmp49 * _tmp91 - _tmp92;
  const Scalar _tmp94 = _tmp65 * _tmp91 - _tmp66 * _tmp93 - _tmp84 * _tmp92 + _tmp87 * _tmp91;
  const Scalar _tmp95 = _tmp68 - _tmp69;
  const Scalar _tmp96 = _tmp72 + _tmp95;
  const Scalar _tmp97 = _tmp71 * _tmp91;
  const Scalar _tmp98 = -_tmp49 * _tmp97 - _tmp75 * _tmp93 + _tmp92 * _tmp96;
  const Scalar _tmp99 = -_tmp81 * _tmp98 - _tmp82 * _tmp93 - _tmp91 * _tmp96 + _tmp97;
  const Scalar _tmp100 = Scalar(1.0) / (_tmp99);
  const Scalar _tmp101 = _tmp100 * _tmp94;
  const Scalar _tmp102 = Scalar(1.0) / (_tmp94);
  const Scalar _tmp103 = _tmp102 * _tmp99;
  const Scalar _tmp104 = _tmp103 * (-_tmp101 * _tmp83 + _tmp49 * _tmp66 - _tmp65);
  const Scalar _tmp105 = _tmp104 + _tmp83;
  const Scalar _tmp106 = _tmp100 * _tmp91;
  const Scalar _tmp107 = _tmp100 * _tmp93;
  const Scalar _tmp108 = _tmp61 * (-_tmp105 * _tmp107 - _tmp49);
  const Scalar _tmp109 = _tmp31 - p_b(1, 0);
  const Scalar _tmp110 = _tmp21 - p_b(0, 0);
  const Scalar _tmp111 =
      std::pow(Scalar(std::pow(_tmp109, Scalar(2)) + std::pow(_tmp110, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp112 = _tmp110 * _tmp111;
  const Scalar _tmp113 = _tmp112 * _tmp64 * (_tmp105 * _tmp106 + _tmp108 * _tmp60 + Scalar(1.0));
  const Scalar _tmp114 = Scalar(1.0) * _tmp75;
  const Scalar _tmp115 = _tmp114 * _tmp81 - Scalar(1.0) * _tmp82;
  const Scalar _tmp116 = _tmp103 * (-_tmp101 * _tmp115 - Scalar(1.0) * _tmp66);
  const Scalar _tmp117 = _tmp115 + _tmp116;
  const Scalar _tmp118 = _tmp61 * (-_tmp107 * _tmp117 + Scalar(1.0));
  const Scalar _tmp119 = _tmp109 * _tmp111;
  const Scalar _tmp120 = _tmp119 * _tmp64 * (_tmp106 * _tmp117 + _tmp118 * _tmp60);
  const Scalar _tmp121 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp122 = _tmp77 * _tmp81 + _tmp80;
  const Scalar _tmp123 = 0;
  const Scalar _tmp124 = _tmp100 * _tmp123;
  const Scalar _tmp125 = _tmp107 * _tmp123 * _tmp61;
  const Scalar _tmp126 = -_tmp112 * _tmp30 + _tmp119 * _tmp20;
  const Scalar _tmp127 = Scalar(1.0) * _tmp102;
  const Scalar _tmp128 = _tmp61 * _tmp93;
  const Scalar _tmp129 = _tmp126 * _tmp64 * (-_tmp127 * _tmp128 * _tmp60 + _tmp127 * _tmp91);
  const Scalar _tmp130 = -_tmp113 * fh1 - _tmp120 * fh1 -
                         _tmp121 * _tmp64 * (_tmp124 * _tmp91 - _tmp125 * _tmp60) - _tmp129 * fh1;
  const Scalar _tmp131 = std::pow(_tmp130, Scalar(-2));
  const Scalar _tmp132 = -_tmp113 - _tmp120 - _tmp129;
  const Scalar _tmp133 = _tmp131 * _tmp132;
  const Scalar _tmp134 = Scalar(0.1034955) * _tmp133;
  const Scalar _tmp135 = _tmp78 + _tmp87;
  const Scalar _tmp136 = _tmp135 * _tmp81;
  const Scalar _tmp137 = Scalar(1.0) / (-_tmp136 + _tmp80 - _tmp84);
  const Scalar _tmp138 = Scalar(1.0) * _tmp137;
  const Scalar _tmp139 = _tmp100 * _tmp98;
  const Scalar _tmp140 = _tmp122 * _tmp137;
  const Scalar _tmp141 = _tmp79 * (-_tmp123 * _tmp139 - _tmp135 * _tmp140 + _tmp78);
  const Scalar _tmp142 = _tmp67 + _tmp95;
  const Scalar _tmp143 = _tmp142 * fh1;
  const Scalar _tmp144 = -_tmp119 * _tmp143 - Scalar(5.1796800000000003) * _tmp25 - _tmp30 * fv1;
  const Scalar _tmp145 = _tmp136 * _tmp138 + Scalar(1.0);
  const Scalar _tmp146 = Scalar(1.0) * _tmp79;
  const Scalar _tmp147 = _tmp138 * _tmp81;
  const Scalar _tmp148 = -Scalar(1.0) * _tmp145 * _tmp146 + Scalar(1.0) * _tmp147;
  const Scalar _tmp149 = _tmp112 * _tmp143 + Scalar(5.1796800000000003) * _tmp17 + _tmp20 * fv1;
  const Scalar _tmp150 = _tmp135 * _tmp79;
  const Scalar _tmp151 = Scalar(1.0) * _tmp138 * _tmp150 - Scalar(1.0) * _tmp138;
  const Scalar _tmp152 = _tmp103 * _tmp138;
  const Scalar _tmp153 = -_tmp127 * _tmp98 + _tmp135 * _tmp152;
  const Scalar _tmp154 = Scalar(1.0) * _tmp126 * (-_tmp146 * _tmp153 + _tmp152);
  const Scalar _tmp155 = _tmp135 * _tmp137;
  const Scalar _tmp156 = -_tmp114 + _tmp116 * _tmp155 - _tmp117 * _tmp139;
  const Scalar _tmp157 = Scalar(1.0) * _tmp119 * (_tmp116 * _tmp138 - _tmp146 * _tmp156);
  const Scalar _tmp158 = _tmp104 * _tmp155 - _tmp105 * _tmp139 + _tmp76;
  const Scalar _tmp159 = Scalar(1.0) * _tmp112 * (_tmp104 * _tmp138 - _tmp146 * _tmp158);
  const Scalar _tmp160 =
      Scalar(1.0) * _tmp121 * (-_tmp122 * _tmp138 - Scalar(1.0) * _tmp141 + Scalar(1.0)) +
      _tmp144 * _tmp148 + _tmp149 * _tmp151 + _tmp154 * fh1 + _tmp157 * fh1 + _tmp159 * fh1;
  const Scalar _tmp161 = Scalar(1.0) / (_tmp130);
  const Scalar _tmp162 = std::asinh(_tmp160 * _tmp161);
  const Scalar _tmp163 = Scalar(1.0) * _tmp162;
  const Scalar _tmp164 = _tmp119 * _tmp142;
  const Scalar _tmp165 = _tmp112 * _tmp142;
  const Scalar _tmp166 = (-_tmp133 * _tmp160 + _tmp161 * (-_tmp148 * _tmp164 + _tmp151 * _tmp165 +
                                                          _tmp154 + _tmp157 + _tmp159)) /
                         std::sqrt(Scalar(_tmp131 * std::pow(_tmp160, Scalar(2)) + 1));
  const Scalar _tmp167 = Scalar(9.6622558468725703) * _tmp162;
  const Scalar _tmp168 =
      -_tmp130 * _tmp167 - std::sqrt(Scalar(std::pow(Scalar(-_tmp43 + p_c(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp46 + p_c(0, 0)), Scalar(2))));
  const Scalar _tmp169 = Scalar(0.1034955) * _tmp161;
  const Scalar _tmp170 = _tmp168 * _tmp169;
  const Scalar _tmp171 = Scalar(9.6622558468725703) * _tmp130;
  const Scalar _tmp172 = _tmp126 * _tmp153 * _tmp79;
  const Scalar _tmp173 = _tmp145 * _tmp79;
  const Scalar _tmp174 = _tmp138 * _tmp149;
  const Scalar _tmp175 = _tmp112 * _tmp79;
  const Scalar _tmp176 = _tmp158 * _tmp175;
  const Scalar _tmp177 = _tmp119 * _tmp156 * _tmp79;
  const Scalar _tmp178 = _tmp121 * _tmp141 + _tmp144 * _tmp173 - _tmp150 * _tmp174 + _tmp172 * fh1 +
                         _tmp176 * fh1 + _tmp177 * fh1;
  const Scalar _tmp179 = _tmp108 * _tmp112;
  const Scalar _tmp180 = _tmp118 * _tmp119;
  const Scalar _tmp181 = _tmp126 * _tmp127;
  const Scalar _tmp182 = _tmp181 * fh1;
  const Scalar _tmp183 = -_tmp121 * _tmp125 - _tmp128 * _tmp182 + _tmp179 * fh1 + _tmp180 * fh1;
  const Scalar _tmp184 = Scalar(1.0) / (_tmp183);
  const Scalar _tmp185 = std::asinh(_tmp178 * _tmp184);
  const Scalar _tmp186 = Scalar(1.0) * _tmp185;
  const Scalar _tmp187 = Scalar(0.1034955) * _tmp184;
  const Scalar _tmp188 = Scalar(9.6622558468725703) * _tmp183;
  const Scalar _tmp189 =
      -_tmp185 * _tmp188 - std::sqrt(Scalar(std::pow(Scalar(-_tmp52 + p_d(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp56 + p_d(0, 0)), Scalar(2))));
  const Scalar _tmp190 = _tmp187 * _tmp189;
  const Scalar _tmp191 = -_tmp128 * _tmp181 + _tmp179 + _tmp180;
  const Scalar _tmp192 = Scalar(9.6622558468725703) * _tmp191;
  const Scalar _tmp193 = std::pow(_tmp183, Scalar(-2));
  const Scalar _tmp194 = _tmp191 * _tmp193;
  const Scalar _tmp195 =
      (-_tmp178 * _tmp194 + _tmp184 * (-_tmp135 * _tmp138 * _tmp142 * _tmp175 - _tmp164 * _tmp173 +
                                       _tmp172 + _tmp176 + _tmp177)) /
      std::sqrt(Scalar(std::pow(_tmp178, Scalar(2)) * _tmp193 + 1));
  const Scalar _tmp196 = Scalar(0.1034955) * _tmp194;
  const Scalar _tmp197 = _tmp100 * _tmp117 * _tmp119;
  const Scalar _tmp198 = _tmp100 * _tmp105 * _tmp112;
  const Scalar _tmp199 = _tmp181 + _tmp197 + _tmp198;
  const Scalar _tmp200 = _tmp121 * _tmp124 + _tmp182 + _tmp197 * fh1 + _tmp198 * fh1;
  const Scalar _tmp201 = Scalar(1.0) / (_tmp200);
  const Scalar _tmp202 = _tmp116 * _tmp119 * _tmp137;
  const Scalar _tmp203 = _tmp126 * _tmp152;
  const Scalar _tmp204 = _tmp104 * _tmp112 * _tmp137;
  const Scalar _tmp205 = _tmp121 * _tmp140 - _tmp144 * _tmp147 + _tmp174 - _tmp202 * fh1 -
                         _tmp203 * fh1 - _tmp204 * fh1;
  const Scalar _tmp206 = std::asinh(_tmp201 * _tmp205);
  const Scalar _tmp207 = Scalar(9.6622558468725703) * _tmp206;
  const Scalar _tmp208 = Scalar(9.6622558468725703) * _tmp200;
  const Scalar _tmp209 = std::pow(_tmp200, Scalar(-2));
  const Scalar _tmp210 = _tmp199 * _tmp209;
  const Scalar _tmp211 =
      (_tmp201 * (_tmp138 * _tmp165 + _tmp147 * _tmp164 - _tmp202 - _tmp203 - _tmp204) -
       _tmp205 * _tmp210) /
      std::sqrt(Scalar(std::pow(_tmp205, Scalar(2)) * _tmp209 + 1));
  const Scalar _tmp212 = Scalar(0.1034955) * _tmp201;
  const Scalar _tmp213 =
      -_tmp200 * _tmp207 - std::sqrt(Scalar(std::pow(Scalar(-_tmp85 + p_a(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp88 + p_a(1, 0)), Scalar(2))));
  const Scalar _tmp214 = Scalar(0.1034955) * _tmp210;
  const Scalar _tmp215 = _tmp212 * _tmp213;
  const Scalar _tmp216 = Scalar(1.0) * _tmp206;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = -Scalar(9.6622558468725703) * _tmp1 * _tmp32 -
               Scalar(9.6622558468725703) * fh1 *
                   (-_tmp0 * _tmp1 - Scalar(1.0) * _tmp0 * _tmp39 * fv1 * std::sinh(_tmp40) -
                    (-Scalar(0.1034955) * _tmp0 * _tmp36 +
                     _tmp37 * (Scalar(9.6622558468725703) * _tmp33 * _tmp39 - _tmp35)) *
                        std::sinh(_tmp38)) +
               Scalar(9.6622558468725703) * std::cosh(_tmp38) -
               Scalar(9.6622558468725703) * std::cosh(_tmp40);
  _res(1, 0) =
      -Scalar(9.6622558468725703) * _tmp132 *
          (_tmp169 * p_c(2, 0) + std::cosh(_tmp163) - std::cosh(_tmp170)) -
      _tmp171 * (-_tmp134 * p_c(2, 0) + Scalar(1.0) * _tmp166 * std::sinh(_tmp163) -
                 (-_tmp134 * _tmp168 + _tmp169 * (-_tmp132 * _tmp167 - _tmp166 * _tmp171)) *
                     std::sinh(_tmp170));
  _res(2, 0) =
      -_tmp188 * (Scalar(1.0) * _tmp195 * std::sinh(_tmp186) - _tmp196 * p_d(2, 0) -
                  (_tmp187 * (-_tmp185 * _tmp192 - _tmp188 * _tmp195) - _tmp189 * _tmp196) *
                      std::sinh(_tmp190)) -
      _tmp192 * (_tmp187 * p_d(2, 0) + std::cosh(_tmp186) - std::cosh(_tmp190));
  _res(3, 0) = -Scalar(9.6622558468725703) * _tmp199 *
                   (_tmp212 * p_a(2, 0) - std::cosh(_tmp215) + std::cosh(_tmp216)) -
               _tmp208 * (Scalar(1.0) * _tmp211 * std::sinh(_tmp216) - _tmp214 * p_a(2, 0) -
                          (_tmp212 * (-_tmp199 * _tmp207 - _tmp208 * _tmp211) - _tmp213 * _tmp214) *
                              std::sinh(_tmp215));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
