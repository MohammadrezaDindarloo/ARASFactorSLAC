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
 * Symbolic function: IK_residual_func_cost1_wrt_fh1_Nl15
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1WrtFh1Nl15(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const Scalar p_init0, const Scalar p_init1, const Scalar p_init2, const Scalar rot_init_x,
    const Scalar rot_init_y, const Scalar rot_init_z, const Scalar rot_init_w,
    const Scalar epsilon) {
  // Total ops: 639

  // Unused inputs
  (void)p_init2;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();

  // Intermediate terms (219)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = _tmp0 * fv1;
  const Scalar _tmp2 = std::asinh(_tmp1);
  const Scalar _tmp3 = Scalar(1.0) * _tmp2;
  const Scalar _tmp4 = Scalar(9.6622558468725703) * _tmp2;
  const Scalar _tmp5 = _DeltaRot[0] * rot_init_z + _DeltaRot[1] * rot_init_w -
                       _DeltaRot[2] * rot_init_x + _DeltaRot[3] * rot_init_y;
  const Scalar _tmp6 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp7 = -_DeltaRot[0] * rot_init_y + _DeltaRot[1] * rot_init_x +
                       _DeltaRot[2] * rot_init_w + _DeltaRot[3] * rot_init_z;
  const Scalar _tmp8 = 1 - 2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp9 = Scalar(0.20999999999999999) * _tmp6 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp10 = _DeltaRot[0] * rot_init_w - _DeltaRot[1] * rot_init_z +
                        _DeltaRot[2] * rot_init_y + _DeltaRot[3] * rot_init_x;
  const Scalar _tmp11 = 2 * _tmp10 * _tmp7;
  const Scalar _tmp12 = -2 * _DeltaRot[0] * rot_init_x - 2 * _DeltaRot[1] * rot_init_y -
                        2 * _DeltaRot[2] * rot_init_z + 2 * _DeltaRot[3] * rot_init_w;
  const Scalar _tmp13 = _tmp12 * _tmp5;
  const Scalar _tmp14 = _tmp11 + _tmp13;
  const Scalar _tmp15 = -Scalar(0.010999999999999999) * _tmp14;
  const Scalar _tmp16 = 2 * _tmp5;
  const Scalar _tmp17 = _tmp10 * _tmp16;
  const Scalar _tmp18 = _tmp12 * _tmp7;
  const Scalar _tmp19 = Scalar(0.20999999999999999) * _tmp17 - Scalar(0.20999999999999999) * _tmp18;
  const Scalar _tmp20 = _tmp15 + _tmp19;
  const Scalar _tmp21 = _tmp20 + _tmp9;
  const Scalar _tmp22 = _tmp21 + p_init0;
  const Scalar _tmp23 = -2 * std::pow(_tmp10, Scalar(2));
  const Scalar _tmp24 = Scalar(0.20999999999999999) * _tmp23 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp25 = _tmp16 * _tmp7;
  const Scalar _tmp26 = _tmp10 * _tmp12;
  const Scalar _tmp27 = _tmp25 - _tmp26;
  const Scalar _tmp28 = -Scalar(0.010999999999999999) * _tmp27;
  const Scalar _tmp29 = Scalar(0.20999999999999999) * _tmp17 + Scalar(0.20999999999999999) * _tmp18;
  const Scalar _tmp30 = _tmp28 + _tmp29;
  const Scalar _tmp31 = _tmp24 + _tmp30;
  const Scalar _tmp32 = _tmp31 + p_init1;
  const Scalar _tmp33 =
      -_tmp4 * fh1 - std::sqrt(Scalar(std::pow(Scalar(-_tmp22 + p_c(0, 0)), Scalar(2)) +
                                      std::pow(Scalar(-_tmp32 + p_c(1, 0)), Scalar(2))));
  const Scalar _tmp34 = Scalar(0.1034955) * _tmp0;
  const Scalar _tmp35 = _tmp33 * _tmp34;
  const Scalar _tmp36 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp37 = Scalar(0.1034955) * _tmp36;
  const Scalar _tmp38 =
      std::pow(Scalar(_tmp36 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp39 = _tmp22 - p_c(0, 0);
  const Scalar _tmp40 = _tmp32 - p_c(1, 0);
  const Scalar _tmp41 = std::pow(Scalar(std::pow(_tmp39, Scalar(2)) + std::pow(_tmp40, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp42 = _tmp39 * _tmp41;
  const Scalar _tmp43 = Scalar(0.20999999999999999) * _tmp11 - Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp44 = -Scalar(0.010999999999999999) * _tmp23 -
                        Scalar(0.010999999999999999) * _tmp6 + Scalar(-0.010999999999999999);
  const Scalar _tmp45 = Scalar(0.20999999999999999) * _tmp25 + Scalar(0.20999999999999999) * _tmp26;
  const Scalar _tmp46 = _tmp44 + _tmp45;
  const Scalar _tmp47 = _tmp43 + _tmp46;
  const Scalar _tmp48 = _tmp47 * fh1;
  const Scalar _tmp49 = Scalar(5.1796800000000003) * _tmp14 + _tmp21 * fv1 + _tmp42 * _tmp48;
  const Scalar _tmp50 = _tmp15 - _tmp19;
  const Scalar _tmp51 = _tmp50 + _tmp9;
  const Scalar _tmp52 = Scalar(1.0) * _tmp51;
  const Scalar _tmp53 = -_tmp9;
  const Scalar _tmp54 = _tmp50 + _tmp53;
  const Scalar _tmp55 = -_tmp24;
  const Scalar _tmp56 = _tmp30 + _tmp55;
  const Scalar _tmp57 = Scalar(1.0) * _tmp56;
  const Scalar _tmp58 = -_tmp57;
  const Scalar _tmp59 = _tmp28 - _tmp29;
  const Scalar _tmp60 = _tmp55 + _tmp59;
  const Scalar _tmp61 = _tmp58 + _tmp60;
  const Scalar _tmp62 = _tmp24 + _tmp59;
  const Scalar _tmp63 = Scalar(1.0) / (_tmp58 + _tmp62);
  const Scalar _tmp64 = _tmp20 + _tmp53;
  const Scalar _tmp65 = _tmp52 - _tmp64;
  const Scalar _tmp66 = _tmp63 * _tmp65;
  const Scalar _tmp67 = _tmp61 * _tmp66;
  const Scalar _tmp68 = Scalar(1.0) / (_tmp52 - _tmp54 - _tmp67);
  const Scalar _tmp69 = Scalar(1.0) * _tmp68;
  const Scalar _tmp70 = _tmp61 * _tmp63;
  const Scalar _tmp71 = Scalar(1.0) * _tmp69 * _tmp70 - Scalar(1.0) * _tmp69;
  const Scalar _tmp72 = _tmp40 * _tmp41;
  const Scalar _tmp73 = -Scalar(5.1796800000000003) * _tmp27 - _tmp31 * fv1 - _tmp48 * _tmp72;
  const Scalar _tmp74 = _tmp67 * _tmp69 + Scalar(1.0);
  const Scalar _tmp75 = Scalar(1.0) * _tmp63;
  const Scalar _tmp76 = _tmp66 * _tmp69;
  const Scalar _tmp77 = -Scalar(1.0) * _tmp74 * _tmp75 + Scalar(1.0) * _tmp76;
  const Scalar _tmp78 = _tmp56 + p_init1;
  const Scalar _tmp79 = _tmp78 - p_b(1, 0);
  const Scalar _tmp80 = _tmp51 + p_init0;
  const Scalar _tmp81 = _tmp80 - p_b(0, 0);
  const Scalar _tmp82 = Scalar(1.0) / (_tmp81);
  const Scalar _tmp83 = _tmp79 * _tmp82;
  const Scalar _tmp84 = -_tmp43;
  const Scalar _tmp85 = _tmp46 + _tmp84;
  const Scalar _tmp86 = _tmp62 + p_init1;
  const Scalar _tmp87 = _tmp86 - p_d(1, 0);
  const Scalar _tmp88 = _tmp64 + p_init0;
  const Scalar _tmp89 = _tmp88 - p_d(0, 0);
  const Scalar _tmp90 = std::pow(Scalar(std::pow(_tmp87, Scalar(2)) + std::pow(_tmp89, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp91 = _tmp89 * _tmp90;
  const Scalar _tmp92 = _tmp44 - _tmp45;
  const Scalar _tmp93 = _tmp43 + _tmp92;
  const Scalar _tmp94 = _tmp87 * _tmp90;
  const Scalar _tmp95 = Scalar(1.0) / (_tmp83 * _tmp91 - _tmp94);
  const Scalar _tmp96 = _tmp95 * (-_tmp85 * _tmp91 + _tmp91 * _tmp93);
  const Scalar _tmp97 = _tmp83 * _tmp93;
  const Scalar _tmp98 = _tmp85 * _tmp94 - _tmp91 * _tmp97;
  const Scalar _tmp99 = _tmp83 * _tmp95;
  const Scalar _tmp100 = _tmp97 + _tmp98 * _tmp99;
  const Scalar _tmp101 = -_tmp100 * _tmp66 + _tmp83 * _tmp96 - _tmp93;
  const Scalar _tmp102 = _tmp54 + p_init0;
  const Scalar _tmp103 = _tmp102 - p_a(0, 0);
  const Scalar _tmp104 = _tmp60 + p_init1;
  const Scalar _tmp105 = _tmp104 - p_a(1, 0);
  const Scalar _tmp106 =
      std::pow(Scalar(std::pow(_tmp103, Scalar(2)) + std::pow(_tmp105, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp107 = _tmp103 * _tmp106;
  const Scalar _tmp108 =
      std::sqrt(Scalar(std::pow(_tmp79, Scalar(2)) + std::pow(_tmp81, Scalar(2))));
  const Scalar _tmp109 = Scalar(1.0) / (_tmp108);
  const Scalar _tmp110 = _tmp108 * _tmp82;
  const Scalar _tmp111 = _tmp110 * (_tmp109 * _tmp51 * _tmp79 - _tmp109 * _tmp56 * _tmp81);
  const Scalar _tmp112 = _tmp105 * _tmp106;
  const Scalar _tmp113 = _tmp111 * _tmp91 + _tmp62 * _tmp91 - _tmp64 * _tmp94;
  const Scalar _tmp114 = _tmp107 * _tmp83 - _tmp112;
  const Scalar _tmp115 = _tmp114 * _tmp95;
  const Scalar _tmp116 =
      _tmp107 * _tmp111 + _tmp107 * _tmp60 - _tmp112 * _tmp54 - _tmp113 * _tmp115;
  const Scalar _tmp117 = _tmp84 + _tmp92;
  const Scalar _tmp118 = -_tmp107 * _tmp97 + _tmp112 * _tmp117 - _tmp115 * _tmp98;
  const Scalar _tmp119 =
      -_tmp107 * _tmp117 + _tmp107 * _tmp93 - _tmp114 * _tmp96 - _tmp118 * _tmp66;
  const Scalar _tmp120 = Scalar(1.0) / (_tmp119);
  const Scalar _tmp121 = _tmp116 * _tmp120;
  const Scalar _tmp122 = Scalar(1.0) / (_tmp116);
  const Scalar _tmp123 = _tmp119 * _tmp122;
  const Scalar _tmp124 = _tmp123 * (-_tmp101 * _tmp121 - _tmp111 + _tmp113 * _tmp99);
  const Scalar _tmp125 = _tmp101 + _tmp124;
  const Scalar _tmp126 = _tmp118 * _tmp120;
  const Scalar _tmp127 = _tmp61 * _tmp68;
  const Scalar _tmp128 = _tmp100 + _tmp124 * _tmp127 - _tmp125 * _tmp126;
  const Scalar _tmp129 = Scalar(1.0) * _tmp42 * (_tmp124 * _tmp69 - _tmp128 * _tmp75);
  const Scalar _tmp130 = Scalar(1.0) * _tmp95;
  const Scalar _tmp131 = _tmp65 * _tmp75 * _tmp95 * _tmp98 - Scalar(1.0) * _tmp96;
  const Scalar _tmp132 = _tmp123 * (-_tmp113 * _tmp130 - _tmp121 * _tmp131);
  const Scalar _tmp133 = _tmp131 + _tmp132;
  const Scalar _tmp134 = -_tmp126 * _tmp133 + _tmp127 * _tmp132 - _tmp130 * _tmp98;
  const Scalar _tmp135 = Scalar(1.0) * _tmp72 * (_tmp132 * _tmp69 - _tmp134 * _tmp75);
  const Scalar _tmp136 = _tmp21 * _tmp72 - _tmp31 * _tmp42;
  const Scalar _tmp137 = Scalar(1.0) * _tmp122;
  const Scalar _tmp138 = _tmp123 * _tmp69;
  const Scalar _tmp139 = -_tmp118 * _tmp137 + _tmp138 * _tmp61;
  const Scalar _tmp140 = Scalar(1.0) * _tmp136 * (_tmp138 - _tmp139 * _tmp75);
  const Scalar _tmp141 = _tmp52 + _tmp57 * _tmp66;
  const Scalar _tmp142 = _tmp141 * _tmp68;
  const Scalar _tmp143 = 0;
  const Scalar _tmp144 = -_tmp126 * _tmp143 - _tmp142 * _tmp61 + _tmp58;
  const Scalar _tmp145 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp146 =
      _tmp129 * fh1 + _tmp135 * fh1 + _tmp140 * fh1 +
      Scalar(1.0) * _tmp145 * (-_tmp141 * _tmp69 - _tmp144 * _tmp75 + Scalar(1.0)) +
      _tmp49 * _tmp71 + _tmp73 * _tmp77;
  const Scalar _tmp147 = _tmp120 * _tmp143;
  const Scalar _tmp148 = _tmp115 * _tmp91;
  const Scalar _tmp149 = _tmp120 * _tmp133;
  const Scalar _tmp150 = _tmp114 * _tmp120;
  const Scalar _tmp151 = _tmp95 * (-_tmp133 * _tmp150 + Scalar(1.0));
  const Scalar _tmp152 = _tmp110 * _tmp72 * (_tmp107 * _tmp149 + _tmp151 * _tmp91);
  const Scalar _tmp153 = _tmp95 * (-_tmp125 * _tmp150 - _tmp83);
  const Scalar _tmp154 = _tmp120 * _tmp125;
  const Scalar _tmp155 = _tmp110 * _tmp42 * (_tmp107 * _tmp154 + _tmp153 * _tmp91 + Scalar(1.0));
  const Scalar _tmp156 = _tmp110 * _tmp136 * (_tmp107 * _tmp137 - _tmp137 * _tmp148);
  const Scalar _tmp157 = -_tmp110 * _tmp145 * (_tmp107 * _tmp147 - _tmp147 * _tmp148) -
                         _tmp152 * fh1 - _tmp155 * fh1 - _tmp156 * fh1;
  const Scalar _tmp158 = Scalar(1.0) / (_tmp157);
  const Scalar _tmp159 = std::asinh(_tmp146 * _tmp158);
  const Scalar _tmp160 = Scalar(9.6622558468725703) * _tmp157;
  const Scalar _tmp161 =
      -_tmp159 * _tmp160 - std::sqrt(Scalar(std::pow(Scalar(-_tmp78 + p_b(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp80 + p_b(0, 0)), Scalar(2))));
  const Scalar _tmp162 = Scalar(0.1034955) * _tmp158;
  const Scalar _tmp163 = _tmp161 * _tmp162;
  const Scalar _tmp164 = Scalar(1.0) * _tmp159;
  const Scalar _tmp165 = -_tmp152 - _tmp155 - _tmp156;
  const Scalar _tmp166 = Scalar(9.6622558468725703) * _tmp165;
  const Scalar _tmp167 = std::pow(_tmp157, Scalar(-2));
  const Scalar _tmp168 = _tmp165 * _tmp167;
  const Scalar _tmp169 = Scalar(0.1034955) * _tmp168;
  const Scalar _tmp170 = _tmp47 * _tmp72;
  const Scalar _tmp171 = _tmp42 * _tmp47;
  const Scalar _tmp172 = (-_tmp146 * _tmp168 + _tmp158 * (_tmp129 + _tmp135 + _tmp140 -
                                                          _tmp170 * _tmp77 + _tmp171 * _tmp71)) /
                         std::sqrt(Scalar(std::pow(_tmp146, Scalar(2)) * _tmp167 + 1));
  const Scalar _tmp173 = _tmp145 * _tmp147;
  const Scalar _tmp174 = _tmp136 * _tmp137;
  const Scalar _tmp175 = _tmp174 * fh1;
  const Scalar _tmp176 = _tmp151 * _tmp72;
  const Scalar _tmp177 = _tmp153 * _tmp42;
  const Scalar _tmp178 = -_tmp115 * _tmp173 - _tmp115 * _tmp175 + _tmp176 * fh1 + _tmp177 * fh1;
  const Scalar _tmp179 = Scalar(1.0) / (_tmp178);
  const Scalar _tmp180 = _tmp63 * _tmp74;
  const Scalar _tmp181 = _tmp49 * _tmp69;
  const Scalar _tmp182 = _tmp42 * _tmp63;
  const Scalar _tmp183 = _tmp128 * _tmp182;
  const Scalar _tmp184 = _tmp134 * _tmp63 * _tmp72;
  const Scalar _tmp185 = _tmp136 * _tmp139 * _tmp63;
  const Scalar _tmp186 = _tmp144 * _tmp145 * _tmp63 + _tmp180 * _tmp73 - _tmp181 * _tmp70 +
                         _tmp183 * fh1 + _tmp184 * fh1 + _tmp185 * fh1;
  const Scalar _tmp187 = std::asinh(_tmp179 * _tmp186);
  const Scalar _tmp188 = Scalar(1.0) * _tmp187;
  const Scalar _tmp189 = Scalar(0.1034955) * _tmp179;
  const Scalar _tmp190 = Scalar(9.6622558468725703) * _tmp178;
  const Scalar _tmp191 =
      -_tmp187 * _tmp190 - std::sqrt(Scalar(std::pow(Scalar(-_tmp86 + p_d(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp88 + p_d(0, 0)), Scalar(2))));
  const Scalar _tmp192 = _tmp189 * _tmp191;
  const Scalar _tmp193 = -_tmp115 * _tmp174 + _tmp176 + _tmp177;
  const Scalar _tmp194 = Scalar(9.6622558468725703) * _tmp193;
  const Scalar _tmp195 = std::pow(_tmp178, Scalar(-2));
  const Scalar _tmp196 = _tmp193 * _tmp195;
  const Scalar _tmp197 = (_tmp179 * (-_tmp170 * _tmp180 - _tmp182 * _tmp47 * _tmp61 * _tmp69 +
                                     _tmp183 + _tmp184 + _tmp185) -
                          _tmp186 * _tmp196) /
                         std::sqrt(Scalar(std::pow(_tmp186, Scalar(2)) * _tmp195 + 1));
  const Scalar _tmp198 = Scalar(0.1034955) * _tmp196;
  const Scalar _tmp199 = _tmp136 * _tmp138;
  const Scalar _tmp200 = _tmp132 * _tmp68 * _tmp72;
  const Scalar _tmp201 = _tmp124 * _tmp42 * _tmp68;
  const Scalar _tmp202 =
      _tmp142 * _tmp145 + _tmp181 - _tmp199 * fh1 - _tmp200 * fh1 - _tmp201 * fh1 - _tmp73 * _tmp76;
  const Scalar _tmp203 = _tmp149 * _tmp72;
  const Scalar _tmp204 = _tmp154 * _tmp42;
  const Scalar _tmp205 = _tmp173 + _tmp175 + _tmp203 * fh1 + _tmp204 * fh1;
  const Scalar _tmp206 = Scalar(1.0) / (_tmp205);
  const Scalar _tmp207 = std::asinh(_tmp202 * _tmp206);
  const Scalar _tmp208 = Scalar(1.0) * _tmp207;
  const Scalar _tmp209 = Scalar(9.6622558468725703) * _tmp205;
  const Scalar _tmp210 =
      -_tmp207 * _tmp209 - std::sqrt(Scalar(std::pow(Scalar(-_tmp102 + p_a(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp104 + p_a(1, 0)), Scalar(2))));
  const Scalar _tmp211 = Scalar(0.1034955) * _tmp206;
  const Scalar _tmp212 = _tmp210 * _tmp211;
  const Scalar _tmp213 = _tmp174 + _tmp203 + _tmp204;
  const Scalar _tmp214 = Scalar(9.6622558468725703) * _tmp213;
  const Scalar _tmp215 = std::pow(_tmp205, Scalar(-2));
  const Scalar _tmp216 = _tmp213 * _tmp215;
  const Scalar _tmp217 = Scalar(0.1034955) * _tmp216;
  const Scalar _tmp218 = (-_tmp202 * _tmp216 + _tmp206 * (_tmp170 * _tmp76 + _tmp171 * _tmp69 -
                                                          _tmp199 - _tmp200 - _tmp201)) /
                         std::sqrt(Scalar(std::pow(_tmp202, Scalar(2)) * _tmp215 + 1));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      -Scalar(9.6622558468725703) * _tmp34 * p_c(2, 0) -
      Scalar(9.6622558468725703) * fh1 *
          (-Scalar(1.0) * _tmp36 * _tmp38 * fv1 * std::sinh(_tmp3) - _tmp37 * p_c(2, 0) -
           (-_tmp33 * _tmp37 + _tmp34 * (Scalar(9.6622558468725703) * _tmp1 * _tmp38 - _tmp4)) *
               std::sinh(_tmp35)) -
      Scalar(9.6622558468725703) * std::cosh(_tmp3) +
      Scalar(9.6622558468725703) * std::cosh(_tmp35);
  _res(1, 0) =
      -_tmp160 * (-_tmp169 * p_b(2, 0) + Scalar(1.0) * _tmp172 * std::sinh(_tmp164) -
                  (-_tmp161 * _tmp169 + _tmp162 * (-_tmp159 * _tmp166 - _tmp160 * _tmp172)) *
                      std::sinh(_tmp163)) -
      _tmp166 * (_tmp162 * p_b(2, 0) - std::cosh(_tmp163) + std::cosh(_tmp164));
  _res(2, 0) =
      -_tmp190 * (Scalar(1.0) * _tmp197 * std::sinh(_tmp188) - _tmp198 * p_d(2, 0) -
                  (_tmp189 * (-_tmp187 * _tmp194 - _tmp190 * _tmp197) - _tmp191 * _tmp198) *
                      std::sinh(_tmp192)) -
      _tmp194 * (_tmp189 * p_d(2, 0) + std::cosh(_tmp188) - std::cosh(_tmp192));
  _res(3, 0) =
      -_tmp209 * (-_tmp217 * p_a(2, 0) + Scalar(1.0) * _tmp218 * std::sinh(_tmp208) -
                  (-_tmp210 * _tmp217 + _tmp211 * (-_tmp207 * _tmp214 - _tmp209 * _tmp218)) *
                      std::sinh(_tmp212)) -
      _tmp214 * (_tmp211 * p_a(2, 0) + std::cosh(_tmp208) - std::cosh(_tmp212));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
