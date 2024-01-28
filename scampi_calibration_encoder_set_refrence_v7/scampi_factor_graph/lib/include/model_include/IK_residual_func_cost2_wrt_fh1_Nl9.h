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
 * Symbolic function: IK_residual_func_cost2_wrt_fh1_Nl9
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2WrtFh1Nl9(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const Scalar p_init0, const Scalar p_init1, const Scalar p_init2, const Scalar rot_init_x,
    const Scalar rot_init_y, const Scalar rot_init_z, const Scalar rot_init_w,
    const Scalar epsilon) {
  // Total ops: 624

  // Unused inputs
  (void)p_init2;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();

  // Intermediate terms (211)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = -_DeltaRot[0] * rot_init_y + _DeltaRot[1] * rot_init_x +
                       _DeltaRot[2] * rot_init_w + _DeltaRot[3] * rot_init_z;
  const Scalar _tmp2 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp3 = _DeltaRot[0] * rot_init_z + _DeltaRot[1] * rot_init_w -
                       _DeltaRot[2] * rot_init_x + _DeltaRot[3] * rot_init_y;
  const Scalar _tmp4 = 1 - 2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp4;
  const Scalar _tmp6 = _DeltaRot[0] * rot_init_w - _DeltaRot[1] * rot_init_z +
                       _DeltaRot[2] * rot_init_y + _DeltaRot[3] * rot_init_x;
  const Scalar _tmp7 = 2 * _tmp1 * _tmp6;
  const Scalar _tmp8 = -2 * _DeltaRot[0] * rot_init_x - 2 * _DeltaRot[1] * rot_init_y -
                       2 * _DeltaRot[2] * rot_init_z + 2 * _DeltaRot[3] * rot_init_w;
  const Scalar _tmp9 = _tmp3 * _tmp8;
  const Scalar _tmp10 = _tmp7 + _tmp9;
  const Scalar _tmp11 = -Scalar(0.010999999999999999) * _tmp10;
  const Scalar _tmp12 = 2 * _tmp3;
  const Scalar _tmp13 = _tmp12 * _tmp6;
  const Scalar _tmp14 = _tmp1 * _tmp8;
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 - Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp11 - _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp5;
  const Scalar _tmp18 = _tmp17 + p_init0;
  const Scalar _tmp19 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp20 = -2 * std::pow(_tmp6, Scalar(2));
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp20 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp22 = _tmp1 * _tmp12;
  const Scalar _tmp23 = _tmp6 * _tmp8;
  const Scalar _tmp24 = _tmp22 - _tmp23;
  const Scalar _tmp25 = -Scalar(0.010999999999999999) * _tmp24;
  const Scalar _tmp26 = -_tmp21 + _tmp25;
  const Scalar _tmp27 = _tmp19 + _tmp26;
  const Scalar _tmp28 = _tmp27 + p_init1;
  const Scalar _tmp29 = _tmp0 * fv1;
  const Scalar _tmp30 = std::asinh(_tmp29);
  const Scalar _tmp31 = Scalar(9.6622558468725703) * _tmp30;
  const Scalar _tmp32 =
      -Scalar(0.1034955) * _tmp31 * fh1 -
      Scalar(0.1034955) * std::sqrt(Scalar(std::pow(Scalar(-_tmp18 + p_b(0, 0)), Scalar(2)) +
                                           std::pow(Scalar(-_tmp28 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp33 = _tmp0 * _tmp32;
  const Scalar _tmp34 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp35 =
      std::pow(Scalar(_tmp34 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp36 = Scalar(1.0) * _tmp30;
  const Scalar _tmp37 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp38 = _tmp11 + _tmp15;
  const Scalar _tmp39 = _tmp38 + _tmp5;
  const Scalar _tmp40 = Scalar(1.0) * _tmp39;
  const Scalar _tmp41 = _tmp21 + _tmp25;
  const Scalar _tmp42 = _tmp19 + _tmp41;
  const Scalar _tmp43 = Scalar(1.0) * _tmp42;
  const Scalar _tmp44 = -_tmp43;
  const Scalar _tmp45 = -_tmp19;
  const Scalar _tmp46 = _tmp41 + _tmp45;
  const Scalar _tmp47 = Scalar(1.0) / (_tmp44 + _tmp46);
  const Scalar _tmp48 = -_tmp5;
  const Scalar _tmp49 = _tmp38 + _tmp48;
  const Scalar _tmp50 = _tmp47 * (_tmp40 - _tmp49);
  const Scalar _tmp51 = _tmp40 + _tmp43 * _tmp50;
  const Scalar _tmp52 = _tmp16 + _tmp48;
  const Scalar _tmp53 = _tmp26 + _tmp45;
  const Scalar _tmp54 = _tmp44 + _tmp53;
  const Scalar _tmp55 = _tmp50 * _tmp54;
  const Scalar _tmp56 = Scalar(1.0) / (_tmp40 - _tmp52 - _tmp55);
  const Scalar _tmp57 = Scalar(1.0) * _tmp56;
  const Scalar _tmp58 = 0;
  const Scalar _tmp59 = Scalar(0.20999999999999999) * _tmp22 + Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp60 = -_tmp59;
  const Scalar _tmp61 =
      -Scalar(0.010999999999999999) * _tmp20 - Scalar(0.010999999999999999) * _tmp4;
  const Scalar _tmp62 = Scalar(0.20999999999999999) * _tmp7 - Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp63 = _tmp61 - _tmp62;
  const Scalar _tmp64 = _tmp60 + _tmp63;
  const Scalar _tmp65 = _tmp52 + p_init0;
  const Scalar _tmp66 = _tmp65 - p_a(0, 0);
  const Scalar _tmp67 = _tmp53 + p_init1;
  const Scalar _tmp68 = _tmp67 - p_a(1, 0);
  const Scalar _tmp69 = std::pow(Scalar(std::pow(_tmp66, Scalar(2)) + std::pow(_tmp68, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp70 = _tmp68 * _tmp69;
  const Scalar _tmp71 = _tmp42 + p_init1;
  const Scalar _tmp72 = _tmp71 - p_c(1, 0);
  const Scalar _tmp73 = _tmp39 + p_init0;
  const Scalar _tmp74 = _tmp73 - p_c(0, 0);
  const Scalar _tmp75 = Scalar(1.0) / (_tmp74);
  const Scalar _tmp76 = _tmp72 * _tmp75;
  const Scalar _tmp77 = _tmp66 * _tmp69;
  const Scalar _tmp78 = -_tmp70 + _tmp76 * _tmp77;
  const Scalar _tmp79 = _tmp59 + _tmp63;
  const Scalar _tmp80 = _tmp46 + p_init1;
  const Scalar _tmp81 = _tmp80 - p_d(1, 0);
  const Scalar _tmp82 = _tmp49 + p_init0;
  const Scalar _tmp83 = _tmp82 - p_d(0, 0);
  const Scalar _tmp84 = std::pow(Scalar(std::pow(_tmp81, Scalar(2)) + std::pow(_tmp83, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp85 = _tmp81 * _tmp84;
  const Scalar _tmp86 = _tmp61 + _tmp62;
  const Scalar _tmp87 = _tmp59 + _tmp86;
  const Scalar _tmp88 = _tmp76 * _tmp87;
  const Scalar _tmp89 = _tmp83 * _tmp84;
  const Scalar _tmp90 = Scalar(1.0) / (_tmp76 * _tmp89 - _tmp85);
  const Scalar _tmp91 = _tmp90 * (_tmp79 * _tmp85 - _tmp88 * _tmp89);
  const Scalar _tmp92 = _tmp77 * _tmp87;
  const Scalar _tmp93 = _tmp64 * _tmp70 - _tmp76 * _tmp92 - _tmp78 * _tmp91;
  const Scalar _tmp94 = _tmp90 * (-_tmp79 * _tmp89 + _tmp87 * _tmp89);
  const Scalar _tmp95 = -_tmp50 * _tmp93 - _tmp64 * _tmp77 - _tmp78 * _tmp94 + _tmp92;
  const Scalar _tmp96 = Scalar(1.0) / (_tmp95);
  const Scalar _tmp97 = _tmp93 * _tmp96;
  const Scalar _tmp98 = _tmp51 * _tmp56;
  const Scalar _tmp99 = _tmp44 - _tmp54 * _tmp98 - _tmp58 * _tmp97;
  const Scalar _tmp100 = Scalar(1.0) * _tmp47;
  const Scalar _tmp101 = _tmp28 - p_b(1, 0);
  const Scalar _tmp102 = _tmp18 - p_b(0, 0);
  const Scalar _tmp103 =
      std::pow(Scalar(std::pow(_tmp101, Scalar(2)) + std::pow(_tmp102, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp104 = _tmp101 * _tmp103;
  const Scalar _tmp105 = _tmp60 + _tmp86;
  const Scalar _tmp106 = _tmp105 * fh1;
  const Scalar _tmp107 = -_tmp104 * _tmp106 - Scalar(5.1796800000000003) * _tmp24 - _tmp27 * fv1;
  const Scalar _tmp108 = _tmp55 * _tmp57 + Scalar(1.0);
  const Scalar _tmp109 = _tmp50 * _tmp57;
  const Scalar _tmp110 = -Scalar(1.0) * _tmp100 * _tmp108 + Scalar(1.0) * _tmp109;
  const Scalar _tmp111 = _tmp102 * _tmp103;
  const Scalar _tmp112 = Scalar(5.1796800000000003) * _tmp10 + _tmp106 * _tmp111 + _tmp17 * fv1;
  const Scalar _tmp113 = _tmp47 * _tmp54;
  const Scalar _tmp114 = Scalar(1.0) * _tmp113 * _tmp57 - Scalar(1.0) * _tmp57;
  const Scalar _tmp115 = _tmp104 * _tmp17 - _tmp111 * _tmp27;
  const Scalar _tmp116 =
      std::sqrt(Scalar(std::pow(_tmp72, Scalar(2)) + std::pow(_tmp74, Scalar(2))));
  const Scalar _tmp117 = Scalar(1.0) / (_tmp116);
  const Scalar _tmp118 = _tmp116 * _tmp75;
  const Scalar _tmp119 = _tmp118 * (_tmp117 * _tmp39 * _tmp72 - _tmp117 * _tmp42 * _tmp74);
  const Scalar _tmp120 = _tmp90 * (_tmp119 * _tmp89 + _tmp46 * _tmp89 - _tmp49 * _tmp85);
  const Scalar _tmp121 = _tmp119 * _tmp77 - _tmp120 * _tmp78 - _tmp52 * _tmp70 + _tmp53 * _tmp77;
  const Scalar _tmp122 = Scalar(1.0) / (_tmp121);
  const Scalar _tmp123 = Scalar(1.0) * _tmp122;
  const Scalar _tmp124 = _tmp122 * _tmp95;
  const Scalar _tmp125 = _tmp124 * _tmp57;
  const Scalar _tmp126 = -_tmp123 * _tmp93 + _tmp125 * _tmp54;
  const Scalar _tmp127 = Scalar(1.0) * _tmp115 * (-_tmp100 * _tmp126 + _tmp125);
  const Scalar _tmp128 = Scalar(1.0) * _tmp91;
  const Scalar _tmp129 = _tmp128 * _tmp50 - Scalar(1.0) * _tmp94;
  const Scalar _tmp130 = _tmp121 * _tmp96;
  const Scalar _tmp131 = _tmp124 * (-Scalar(1.0) * _tmp120 - _tmp129 * _tmp130);
  const Scalar _tmp132 = _tmp54 * _tmp56;
  const Scalar _tmp133 = _tmp129 + _tmp131;
  const Scalar _tmp134 = -_tmp128 + _tmp131 * _tmp132 - _tmp133 * _tmp97;
  const Scalar _tmp135 = Scalar(1.0) * _tmp104 * (-_tmp100 * _tmp134 + _tmp131 * _tmp57);
  const Scalar _tmp136 = _tmp76 * _tmp91 + _tmp88;
  const Scalar _tmp137 = -_tmp136 * _tmp50 + _tmp76 * _tmp94 - _tmp87;
  const Scalar _tmp138 = _tmp124 * (-_tmp119 + _tmp120 * _tmp76 - _tmp130 * _tmp137);
  const Scalar _tmp139 = _tmp96 * (_tmp137 + _tmp138);
  const Scalar _tmp140 = _tmp132 * _tmp138 + _tmp136 - _tmp139 * _tmp93;
  const Scalar _tmp141 = Scalar(1.0) * _tmp111 * (-_tmp100 * _tmp140 + _tmp138 * _tmp57);
  const Scalar _tmp142 = _tmp107 * _tmp110 + _tmp112 * _tmp114 + _tmp127 * fh1 + _tmp135 * fh1 +
                         _tmp141 * fh1 +
                         Scalar(1.0) * _tmp37 * (-_tmp100 * _tmp99 - _tmp51 * _tmp57 + Scalar(1.0));
  const Scalar _tmp143 = _tmp90 * (-_tmp139 * _tmp78 - _tmp76);
  const Scalar _tmp144 = _tmp111 * _tmp118 * (_tmp139 * _tmp77 + _tmp143 * _tmp89 + Scalar(1.0));
  const Scalar _tmp145 = _tmp133 * _tmp96;
  const Scalar _tmp146 = _tmp90 * (-_tmp145 * _tmp78 + Scalar(1.0));
  const Scalar _tmp147 = _tmp104 * _tmp118 * (_tmp145 * _tmp77 + _tmp146 * _tmp89);
  const Scalar _tmp148 = _tmp58 * _tmp96;
  const Scalar _tmp149 = _tmp78 * _tmp90;
  const Scalar _tmp150 = _tmp149 * _tmp89;
  const Scalar _tmp151 = _tmp115 * _tmp118 * (-_tmp123 * _tmp150 + _tmp123 * _tmp77);
  const Scalar _tmp152 = -_tmp118 * _tmp37 * (-_tmp148 * _tmp150 + _tmp148 * _tmp77) -
                         _tmp144 * fh1 - _tmp147 * fh1 - _tmp151 * fh1;
  const Scalar _tmp153 = Scalar(1.0) / (_tmp152);
  const Scalar _tmp154 = std::asinh(_tmp142 * _tmp153);
  const Scalar _tmp155 = Scalar(1.0) * _tmp154;
  const Scalar _tmp156 = std::pow(_tmp152, Scalar(-2));
  const Scalar _tmp157 = _tmp104 * _tmp105;
  const Scalar _tmp158 = _tmp105 * _tmp111;
  const Scalar _tmp159 = -_tmp144 - _tmp147 - _tmp151;
  const Scalar _tmp160 = _tmp156 * _tmp159;
  const Scalar _tmp161 = (-_tmp142 * _tmp160 + _tmp153 * (-_tmp110 * _tmp157 + _tmp114 * _tmp158 +
                                                          _tmp127 + _tmp135 + _tmp141)) /
                         std::sqrt(Scalar(std::pow(_tmp142, Scalar(2)) * _tmp156 + 1));
  const Scalar _tmp162 = Scalar(9.6622558468725703) * _tmp154;
  const Scalar _tmp163 =
      -_tmp152 * _tmp162 - std::sqrt(Scalar(std::pow(Scalar(-_tmp71 + p_c(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp73 + p_c(0, 0)), Scalar(2))));
  const Scalar _tmp164 = Scalar(0.1034955) * _tmp153;
  const Scalar _tmp165 = _tmp163 * _tmp164;
  const Scalar _tmp166 = Scalar(9.6622558468725703) * _tmp152;
  const Scalar _tmp167 = _tmp115 * _tmp126 * _tmp47;
  const Scalar _tmp168 = _tmp108 * _tmp47;
  const Scalar _tmp169 = _tmp112 * _tmp57;
  const Scalar _tmp170 = _tmp111 * _tmp47;
  const Scalar _tmp171 = _tmp140 * _tmp170;
  const Scalar _tmp172 = _tmp104 * _tmp134 * _tmp47;
  const Scalar _tmp173 = _tmp107 * _tmp168 - _tmp113 * _tmp169 + _tmp167 * fh1 + _tmp171 * fh1 +
                         _tmp172 * fh1 + _tmp37 * _tmp47 * _tmp99;
  const Scalar _tmp174 = _tmp111 * _tmp143;
  const Scalar _tmp175 = _tmp104 * _tmp146;
  const Scalar _tmp176 = _tmp115 * _tmp123;
  const Scalar _tmp177 = _tmp176 * fh1;
  const Scalar _tmp178 = _tmp148 * _tmp37;
  const Scalar _tmp179 = -_tmp149 * _tmp177 - _tmp149 * _tmp178 + _tmp174 * fh1 + _tmp175 * fh1;
  const Scalar _tmp180 = Scalar(1.0) / (_tmp179);
  const Scalar _tmp181 = std::asinh(_tmp173 * _tmp180);
  const Scalar _tmp182 = Scalar(1.0) * _tmp181;
  const Scalar _tmp183 = std::pow(_tmp179, Scalar(-2));
  const Scalar _tmp184 = -_tmp149 * _tmp176 + _tmp174 + _tmp175;
  const Scalar _tmp185 = _tmp183 * _tmp184;
  const Scalar _tmp186 =
      (-_tmp173 * _tmp185 + _tmp180 * (-_tmp105 * _tmp170 * _tmp54 * _tmp57 - _tmp157 * _tmp168 +
                                       _tmp167 + _tmp171 + _tmp172)) /
      std::sqrt(Scalar(std::pow(_tmp173, Scalar(2)) * _tmp183 + 1));
  const Scalar _tmp187 = Scalar(9.6622558468725703) * _tmp179;
  const Scalar _tmp188 = Scalar(9.6622558468725703) * _tmp184;
  const Scalar _tmp189 = Scalar(0.1034955) * _tmp180;
  const Scalar _tmp190 =
      -_tmp181 * _tmp187 - std::sqrt(Scalar(std::pow(Scalar(-_tmp80 + p_d(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp82 + p_d(0, 0)), Scalar(2))));
  const Scalar _tmp191 = _tmp189 * _tmp190;
  const Scalar _tmp192 = _tmp104 * _tmp145;
  const Scalar _tmp193 = _tmp111 * _tmp139;
  const Scalar _tmp194 = _tmp176 + _tmp192 + _tmp193;
  const Scalar _tmp195 = _tmp177 + _tmp178 + _tmp192 * fh1 + _tmp193 * fh1;
  const Scalar _tmp196 = Scalar(1.0) / (_tmp195);
  const Scalar _tmp197 = _tmp104 * _tmp131 * _tmp56;
  const Scalar _tmp198 = _tmp115 * _tmp125;
  const Scalar _tmp199 = _tmp111 * _tmp138 * _tmp56;
  const Scalar _tmp200 = -_tmp107 * _tmp109 + _tmp169 - _tmp197 * fh1 - _tmp198 * fh1 -
                         _tmp199 * fh1 + _tmp37 * _tmp98;
  const Scalar _tmp201 = std::asinh(_tmp196 * _tmp200);
  const Scalar _tmp202 = Scalar(9.6622558468725703) * _tmp201;
  const Scalar _tmp203 = Scalar(9.6622558468725703) * _tmp195;
  const Scalar _tmp204 = std::pow(_tmp195, Scalar(-2));
  const Scalar _tmp205 = _tmp194 * _tmp204;
  const Scalar _tmp206 =
      (_tmp196 * (_tmp109 * _tmp157 + _tmp158 * _tmp57 - _tmp197 - _tmp198 - _tmp199) -
       _tmp200 * _tmp205) /
      std::sqrt(Scalar(std::pow(_tmp200, Scalar(2)) * _tmp204 + 1));
  const Scalar _tmp207 = Scalar(0.1034955) * _tmp196;
  const Scalar _tmp208 =
      -_tmp195 * _tmp202 - std::sqrt(Scalar(std::pow(Scalar(-_tmp65 + p_a(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp67 + p_a(1, 0)), Scalar(2))));
  const Scalar _tmp209 = _tmp207 * _tmp208;
  const Scalar _tmp210 = Scalar(1.0) * _tmp201;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      Scalar(9.6622558468725703) * fh1 *
          (Scalar(1.0) * _tmp34 * _tmp35 * fv1 * std::cosh(_tmp36) -
           (Scalar(0.1034955) * _tmp0 * (Scalar(9.6622558468725703) * _tmp29 * _tmp35 - _tmp31) -
            _tmp32 * _tmp34) *
               std::cosh(_tmp33)) -
      Scalar(9.6622558468725703) * std::sinh(_tmp33) -
      Scalar(9.6622558468725703) * std::sinh(_tmp36);
  _res(1, 0) = Scalar(9.6622558468725703) * _tmp159 * (-std::sinh(_tmp155) - std::sinh(_tmp165)) +
               _tmp166 * (-Scalar(1.0) * _tmp161 * std::cosh(_tmp155) -
                          (-Scalar(0.1034955) * _tmp160 * _tmp163 +
                           _tmp164 * (-_tmp159 * _tmp162 - _tmp161 * _tmp166)) *
                              std::cosh(_tmp165));
  _res(2, 0) = _tmp187 * (-Scalar(1.0) * _tmp186 * std::cosh(_tmp182) -
                          (-Scalar(0.1034955) * _tmp185 * _tmp190 +
                           _tmp189 * (-_tmp181 * _tmp188 - _tmp186 * _tmp187)) *
                              std::cosh(_tmp191)) +
               _tmp188 * (-std::sinh(_tmp182) - std::sinh(_tmp191));
  _res(3, 0) = Scalar(9.6622558468725703) * _tmp194 * (-std::sinh(_tmp209) - std::sinh(_tmp210)) +
               _tmp203 * (-Scalar(1.0) * _tmp206 * std::cosh(_tmp210) -
                          (-Scalar(0.1034955) * _tmp205 * _tmp208 +
                           _tmp207 * (-_tmp194 * _tmp202 - _tmp203 * _tmp206)) *
                              std::cosh(_tmp209));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
