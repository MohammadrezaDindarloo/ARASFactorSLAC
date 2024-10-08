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
 * Symbolic function: IK_residual_func_cost2_wrt_fh1_Nl19
 *
 * Args:
 *     fh1: Scalar
 *     fv1: Scalar
 *     DeltaRot: Rot3
 *     position_vector: Matrix31
 *     encoder: Matrix41
 *     p_a: Matrix31
 *     p_b: Matrix31
 *     p_c: Matrix31
 *     p_d: Matrix31
 *     Rot_init: Rot3
 *     epsilon: Scalar
 *
 * Outputs:
 *     res: Matrix41
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2WrtFh1Nl19(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 620

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (211)
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
  const Scalar _tmp11 = -_tmp10 + _tmp9;
  const Scalar _tmp12 = -Scalar(0.010999999999999999) * _tmp11;
  const Scalar _tmp13 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp14 = 1 - 2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp12 + _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp8;
  const Scalar _tmp18 = _tmp17 + position_vector(1, 0);
  const Scalar _tmp19 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp20 = Scalar(0.20999999999999999) * _tmp14 + Scalar(0.20999999999999999) * _tmp19;
  const Scalar _tmp21 = -_tmp20;
  const Scalar _tmp22 = 2 * _tmp0 * _tmp4;
  const Scalar _tmp23 = _tmp1 * _tmp5;
  const Scalar _tmp24 = _tmp22 + _tmp23;
  const Scalar _tmp25 = -Scalar(0.010999999999999999) * _tmp24;
  const Scalar _tmp26 = Scalar(0.20999999999999999) * _tmp3 - Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp27 = _tmp25 + _tmp26;
  const Scalar _tmp28 = _tmp21 + _tmp27;
  const Scalar _tmp29 = _tmp28 + position_vector(0, 0);
  const Scalar _tmp30 = Scalar(1.0) / (fh1);
  const Scalar _tmp31 = _tmp30 * fv1;
  const Scalar _tmp32 = std::asinh(_tmp31);
  const Scalar _tmp33 = Scalar(1.4083112389913199) * _tmp32;
  const Scalar _tmp34 =
      -_tmp33 * fh1 - std::sqrt(Scalar(std::pow(Scalar(-_tmp18 + p_d(1, 0)), Scalar(2)) +
                                       std::pow(Scalar(-_tmp29 + p_d(0, 0)), Scalar(2))));
  const Scalar _tmp35 = Scalar(0.71007031138673404) * _tmp30;
  const Scalar _tmp36 = _tmp34 * _tmp35;
  const Scalar _tmp37 = Scalar(1.0) * _tmp32;
  const Scalar _tmp38 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp39 =
      std::pow(Scalar(_tmp38 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp40 = Scalar(0.20999999999999999) * _tmp22 - Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp41 = -Scalar(0.010999999999999999) * _tmp13 -
                        Scalar(0.010999999999999999) * _tmp19 + Scalar(-0.010999999999999999);
  const Scalar _tmp42 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp43 = _tmp41 - _tmp42;
  const Scalar _tmp44 = _tmp40 + _tmp43;
  const Scalar _tmp45 = _tmp25 - _tmp26;
  const Scalar _tmp46 = _tmp20 + _tmp45;
  const Scalar _tmp47 = _tmp46 + position_vector(0, 0);
  const Scalar _tmp48 = _tmp47 - p_b(0, 0);
  const Scalar _tmp49 = _tmp12 - _tmp15;
  const Scalar _tmp50 = _tmp49 + _tmp7;
  const Scalar _tmp51 = _tmp50 + position_vector(1, 0);
  const Scalar _tmp52 = _tmp51 - p_b(1, 0);
  const Scalar _tmp53 = std::pow(Scalar(std::pow(_tmp48, Scalar(2)) + std::pow(_tmp52, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp54 = _tmp48 * _tmp53;
  const Scalar _tmp55 = _tmp52 * _tmp53;
  const Scalar _tmp56 = _tmp41 + _tmp42;
  const Scalar _tmp57 = _tmp40 + _tmp56;
  const Scalar _tmp58 = _tmp16 + _tmp7;
  const Scalar _tmp59 = _tmp58 + position_vector(1, 0);
  const Scalar _tmp60 = _tmp59 - p_c(1, 0);
  const Scalar _tmp61 = _tmp20 + _tmp27;
  const Scalar _tmp62 = _tmp61 + position_vector(0, 0);
  const Scalar _tmp63 = _tmp62 - p_c(0, 0);
  const Scalar _tmp64 = std::pow(Scalar(std::pow(_tmp60, Scalar(2)) + std::pow(_tmp63, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp65 = _tmp60 * _tmp64;
  const Scalar _tmp66 = -_tmp40;
  const Scalar _tmp67 = _tmp43 + _tmp66;
  const Scalar _tmp68 = _tmp21 + _tmp45;
  const Scalar _tmp69 = _tmp68 + position_vector(0, 0);
  const Scalar _tmp70 = _tmp69 - p_a(0, 0);
  const Scalar _tmp71 = Scalar(1.0) / (_tmp70);
  const Scalar _tmp72 = _tmp49 + _tmp8;
  const Scalar _tmp73 = _tmp72 + position_vector(1, 0);
  const Scalar _tmp74 = _tmp73 - p_a(1, 0);
  const Scalar _tmp75 = _tmp71 * _tmp74;
  const Scalar _tmp76 = _tmp67 * _tmp75;
  const Scalar _tmp77 = _tmp63 * _tmp64;
  const Scalar _tmp78 = _tmp57 * _tmp65 - _tmp76 * _tmp77;
  const Scalar _tmp79 = _tmp54 * _tmp75 - _tmp55;
  const Scalar _tmp80 = Scalar(1.0) / (-_tmp65 + _tmp75 * _tmp77);
  const Scalar _tmp81 = _tmp79 * _tmp80;
  const Scalar _tmp82 = _tmp54 * _tmp67;
  const Scalar _tmp83 = _tmp44 * _tmp55 - _tmp75 * _tmp82 - _tmp78 * _tmp81;
  const Scalar _tmp84 = Scalar(1.0) * _tmp72;
  const Scalar _tmp85 = -_tmp84;
  const Scalar _tmp86 = Scalar(1.0) / (_tmp58 + _tmp85);
  const Scalar _tmp87 = Scalar(1.0) * _tmp68;
  const Scalar _tmp88 = -_tmp61 + _tmp87;
  const Scalar _tmp89 = _tmp86 * _tmp88;
  const Scalar _tmp90 = -_tmp57 * _tmp77 + _tmp67 * _tmp77;
  const Scalar _tmp91 = -_tmp44 * _tmp54 - _tmp81 * _tmp90 + _tmp82 - _tmp83 * _tmp89;
  const Scalar _tmp92 = Scalar(1.0) / (_tmp91);
  const Scalar _tmp93 =
      std::sqrt(Scalar(std::pow(_tmp70, Scalar(2)) + std::pow(_tmp74, Scalar(2))));
  const Scalar _tmp94 = Scalar(1.0) / (_tmp93);
  const Scalar _tmp95 = _tmp71 * _tmp93;
  const Scalar _tmp96 = _tmp95 * (_tmp68 * _tmp74 * _tmp94 - _tmp70 * _tmp72 * _tmp94);
  const Scalar _tmp97 = _tmp58 * _tmp77 - _tmp61 * _tmp65 + _tmp77 * _tmp96;
  const Scalar _tmp98 = _tmp75 * _tmp80;
  const Scalar _tmp99 = _tmp76 + _tmp78 * _tmp98;
  const Scalar _tmp100 = -_tmp67 - _tmp89 * _tmp99 + _tmp90 * _tmp98;
  const Scalar _tmp101 = -_tmp46 * _tmp55 + _tmp50 * _tmp54 + _tmp54 * _tmp96 - _tmp81 * _tmp97;
  const Scalar _tmp102 = _tmp101 * _tmp92;
  const Scalar _tmp103 = Scalar(1.0) / (_tmp101);
  const Scalar _tmp104 = _tmp103 * _tmp91;
  const Scalar _tmp105 = _tmp104 * (-_tmp100 * _tmp102 - _tmp96 + _tmp97 * _tmp98);
  const Scalar _tmp106 = _tmp92 * (_tmp100 + _tmp105);
  const Scalar _tmp107 = _tmp80 * (-_tmp106 * _tmp79 - _tmp75);
  const Scalar _tmp108 = _tmp29 - p_d(0, 0);
  const Scalar _tmp109 = _tmp18 - p_d(1, 0);
  const Scalar _tmp110 =
      std::pow(Scalar(std::pow(_tmp108, Scalar(2)) + std::pow(_tmp109, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp111 = _tmp108 * _tmp110;
  const Scalar _tmp112 = _tmp111 * _tmp95 * (_tmp106 * _tmp54 + _tmp107 * _tmp77 + Scalar(1.0));
  const Scalar _tmp113 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp114 = _tmp84 * _tmp89 + _tmp87;
  const Scalar _tmp115 = 0;
  const Scalar _tmp116 = _tmp77 * _tmp81;
  const Scalar _tmp117 = _tmp109 * _tmp110;
  const Scalar _tmp118 = -_tmp111 * _tmp17 + _tmp117 * _tmp28;
  const Scalar _tmp119 = Scalar(1.0) * _tmp103;
  const Scalar _tmp120 = _tmp118 * _tmp95 * (-_tmp116 * _tmp119 + _tmp119 * _tmp54);
  const Scalar _tmp121 = Scalar(1.0) * _tmp80;
  const Scalar _tmp122 = Scalar(1.0) * _tmp86;
  const Scalar _tmp123 = -_tmp121 * _tmp90 + _tmp122 * _tmp78 * _tmp80 * _tmp88;
  const Scalar _tmp124 = _tmp104 * (-_tmp102 * _tmp123 - _tmp121 * _tmp97);
  const Scalar _tmp125 = _tmp92 * (_tmp123 + _tmp124);
  const Scalar _tmp126 = _tmp80 * (-_tmp125 * _tmp79 + Scalar(1.0));
  const Scalar _tmp127 = _tmp117 * _tmp95 * (_tmp125 * _tmp54 + _tmp126 * _tmp77);
  const Scalar _tmp128 = -_tmp112 * fh1 -
                         _tmp113 * _tmp95 * (-_tmp115 * _tmp116 + _tmp115 * _tmp54) -
                         _tmp120 * fh1 - _tmp127 * fh1;
  const Scalar _tmp129 = Scalar(1.0) / (_tmp128);
  const Scalar _tmp130 = _tmp56 + _tmp66;
  const Scalar _tmp131 = _tmp130 * fh1;
  const Scalar _tmp132 = -Scalar(40.024799999999999) * _tmp11 - _tmp117 * _tmp131 - _tmp17 * fv1;
  const Scalar _tmp133 = _tmp50 + _tmp85;
  const Scalar _tmp134 = _tmp133 * _tmp89;
  const Scalar _tmp135 = Scalar(1.0) / (-_tmp134 - _tmp46 + _tmp87);
  const Scalar _tmp136 = Scalar(1.0) * _tmp135;
  const Scalar _tmp137 = _tmp134 * _tmp136 + Scalar(1.0);
  const Scalar _tmp138 = _tmp136 * _tmp89;
  const Scalar _tmp139 = -Scalar(1.0) * _tmp122 * _tmp137 + Scalar(1.0) * _tmp138;
  const Scalar _tmp140 = _tmp104 * _tmp136;
  const Scalar _tmp141 = -_tmp119 * _tmp83 + _tmp133 * _tmp140;
  const Scalar _tmp142 = Scalar(1.0) * _tmp118 * (-_tmp122 * _tmp141 + _tmp140);
  const Scalar _tmp143 = _tmp133 * _tmp135;
  const Scalar _tmp144 = -_tmp121 * _tmp78 + _tmp124 * _tmp143 - _tmp125 * _tmp83;
  const Scalar _tmp145 = Scalar(1.0) * _tmp117 * (-_tmp122 * _tmp144 + _tmp124 * _tmp136);
  const Scalar _tmp146 = _tmp105 * _tmp143 - _tmp106 * _tmp83 + _tmp99;
  const Scalar _tmp147 = Scalar(1.0) * _tmp111 * (_tmp105 * _tmp136 - _tmp122 * _tmp146);
  const Scalar _tmp148 = _tmp111 * _tmp131 + Scalar(40.024799999999999) * _tmp24 + _tmp28 * fv1;
  const Scalar _tmp149 = _tmp133 * _tmp86;
  const Scalar _tmp150 = Scalar(1.0) * _tmp136 * _tmp149 - Scalar(1.0) * _tmp136;
  const Scalar _tmp151 = _tmp114 * _tmp135;
  const Scalar _tmp152 = -_tmp115 * _tmp83 - _tmp133 * _tmp151 + _tmp85;
  const Scalar _tmp153 =
      Scalar(1.0) * _tmp113 * (-_tmp114 * _tmp136 - _tmp122 * _tmp152 + Scalar(1.0)) +
      _tmp132 * _tmp139 + _tmp142 * fh1 + _tmp145 * fh1 + _tmp147 * fh1 + _tmp148 * _tmp150;
  const Scalar _tmp154 = std::asinh(_tmp129 * _tmp153);
  const Scalar _tmp155 = Scalar(1.4083112389913199) * _tmp128;
  const Scalar _tmp156 =
      -_tmp154 * _tmp155 - std::sqrt(Scalar(std::pow(Scalar(-_tmp69 + p_a(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp73 + p_a(1, 0)), Scalar(2))));
  const Scalar _tmp157 = Scalar(0.71007031138673404) * _tmp129;
  const Scalar _tmp158 = _tmp156 * _tmp157;
  const Scalar _tmp159 = Scalar(1.0) * _tmp154;
  const Scalar _tmp160 = -_tmp112 - _tmp120 - _tmp127;
  const Scalar _tmp161 = Scalar(1.4083112389913199) * _tmp160;
  const Scalar _tmp162 = _tmp117 * _tmp130;
  const Scalar _tmp163 = _tmp111 * _tmp130;
  const Scalar _tmp164 = std::pow(_tmp128, Scalar(-2));
  const Scalar _tmp165 = _tmp160 * _tmp164;
  const Scalar _tmp166 =
      (_tmp129 * (-_tmp139 * _tmp162 + _tmp142 + _tmp145 + _tmp147 + _tmp150 * _tmp163) -
       _tmp153 * _tmp165) /
      std::sqrt(Scalar(std::pow(_tmp153, Scalar(2)) * _tmp164 + 1));
  const Scalar _tmp167 = _tmp117 * _tmp126;
  const Scalar _tmp168 = _tmp118 * _tmp119;
  const Scalar _tmp169 = _tmp168 * fh1;
  const Scalar _tmp170 = _tmp107 * _tmp111;
  const Scalar _tmp171 = _tmp113 * _tmp115;
  const Scalar _tmp172 = _tmp167 * fh1 - _tmp169 * _tmp81 + _tmp170 * fh1 - _tmp171 * _tmp81;
  const Scalar _tmp173 = Scalar(1.0) / (_tmp172);
  const Scalar _tmp174 = _tmp118 * _tmp141 * _tmp86;
  const Scalar _tmp175 = _tmp137 * _tmp86;
  const Scalar _tmp176 = _tmp117 * _tmp144 * _tmp86;
  const Scalar _tmp177 = _tmp136 * _tmp148;
  const Scalar _tmp178 = _tmp111 * _tmp86;
  const Scalar _tmp179 = _tmp146 * _tmp178;
  const Scalar _tmp180 = _tmp113 * _tmp152 * _tmp86 + _tmp132 * _tmp175 - _tmp149 * _tmp177 +
                         _tmp174 * fh1 + _tmp176 * fh1 + _tmp179 * fh1;
  const Scalar _tmp181 = std::asinh(_tmp173 * _tmp180);
  const Scalar _tmp182 = Scalar(1.4083112389913199) * _tmp172;
  const Scalar _tmp183 =
      -_tmp181 * _tmp182 - std::sqrt(Scalar(std::pow(Scalar(-_tmp59 + p_c(1, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp62 + p_c(0, 0)), Scalar(2))));
  const Scalar _tmp184 = Scalar(0.71007031138673404) * _tmp173;
  const Scalar _tmp185 = _tmp183 * _tmp184;
  const Scalar _tmp186 = Scalar(1.0) * _tmp181;
  const Scalar _tmp187 = _tmp167 - _tmp168 * _tmp81 + _tmp170;
  const Scalar _tmp188 = Scalar(1.4083112389913199) * _tmp187;
  const Scalar _tmp189 = std::pow(_tmp172, Scalar(-2));
  const Scalar _tmp190 = _tmp187 * _tmp189;
  const Scalar _tmp191 = (_tmp173 * (-_tmp130 * _tmp133 * _tmp136 * _tmp178 - _tmp162 * _tmp175 +
                                     _tmp174 + _tmp176 + _tmp179) -
                          _tmp180 * _tmp190) /
                         std::sqrt(Scalar(std::pow(_tmp180, Scalar(2)) * _tmp189 + 1));
  const Scalar _tmp192 = _tmp117 * _tmp125;
  const Scalar _tmp193 = _tmp106 * _tmp111;
  const Scalar _tmp194 = _tmp169 + _tmp171 + _tmp192 * fh1 + _tmp193 * fh1;
  const Scalar _tmp195 = Scalar(1.0) / (_tmp194);
  const Scalar _tmp196 = _tmp105 * _tmp111 * _tmp135;
  const Scalar _tmp197 = _tmp118 * _tmp140;
  const Scalar _tmp198 = _tmp117 * _tmp124 * _tmp135;
  const Scalar _tmp199 = _tmp113 * _tmp151 - _tmp132 * _tmp138 + _tmp177 - _tmp196 * fh1 -
                         _tmp197 * fh1 - _tmp198 * fh1;
  const Scalar _tmp200 = std::asinh(_tmp195 * _tmp199);
  const Scalar _tmp201 = Scalar(1.4083112389913199) * _tmp194;
  const Scalar _tmp202 =
      -_tmp200 * _tmp201 - std::sqrt(Scalar(std::pow(Scalar(-_tmp47 + p_b(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp51 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp203 = std::pow(_tmp194, Scalar(-2));
  const Scalar _tmp204 = _tmp168 + _tmp192 + _tmp193;
  const Scalar _tmp205 = _tmp203 * _tmp204;
  const Scalar _tmp206 = Scalar(1.4083112389913199) * _tmp204;
  const Scalar _tmp207 =
      (_tmp195 * (_tmp136 * _tmp163 + _tmp138 * _tmp162 - _tmp196 - _tmp197 - _tmp198) -
       _tmp199 * _tmp205) /
      std::sqrt(Scalar(std::pow(_tmp199, Scalar(2)) * _tmp203 + 1));
  const Scalar _tmp208 = Scalar(0.71007031138673404) * _tmp195;
  const Scalar _tmp209 = _tmp202 * _tmp208;
  const Scalar _tmp210 = Scalar(1.0) * _tmp200;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = Scalar(1.4083112389913199) * fh1 *
                   (Scalar(1.0) * _tmp38 * _tmp39 * fv1 * std::cosh(_tmp37) -
                    (-Scalar(0.71007031138673404) * _tmp34 * _tmp38 +
                     _tmp35 * (Scalar(1.4083112389913199) * _tmp31 * _tmp39 - _tmp33)) *
                        std::cosh(_tmp36)) -
               Scalar(1.4083112389913199) * std::sinh(_tmp36) -
               Scalar(1.4083112389913199) * std::sinh(_tmp37);
  _res(1, 0) = _tmp155 * (-Scalar(1.0) * _tmp166 * std::cosh(_tmp159) -
                          (-Scalar(0.71007031138673404) * _tmp156 * _tmp165 +
                           _tmp157 * (-_tmp154 * _tmp161 - _tmp155 * _tmp166)) *
                              std::cosh(_tmp158)) +
               _tmp161 * (-std::sinh(_tmp158) - std::sinh(_tmp159));
  _res(2, 0) = _tmp182 * (-Scalar(1.0) * _tmp191 * std::cosh(_tmp186) -
                          (-Scalar(0.71007031138673404) * _tmp183 * _tmp190 +
                           _tmp184 * (-_tmp181 * _tmp188 - _tmp182 * _tmp191)) *
                              std::cosh(_tmp185)) +
               _tmp188 * (-std::sinh(_tmp185) - std::sinh(_tmp186));
  _res(3, 0) = _tmp201 * (-Scalar(1.0) * _tmp207 * std::cosh(_tmp210) -
                          (-Scalar(0.71007031138673404) * _tmp202 * _tmp205 +
                           _tmp208 * (-_tmp200 * _tmp206 - _tmp201 * _tmp207)) *
                              std::cosh(_tmp209)) +
               _tmp206 * (-std::sinh(_tmp209) - std::sinh(_tmp210));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
