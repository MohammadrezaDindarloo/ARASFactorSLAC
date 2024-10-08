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
 * Symbolic function: IK_residual_func_cost2_wrt_fh1_Nl12
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2WrtFh1Nl12(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 618

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (211)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp1 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp2 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp3 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp5 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp6 = 2 * _tmp2 * _tmp5;
  const Scalar _tmp7 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp8 = _tmp0 * _tmp7;
  const Scalar _tmp9 = _tmp6 - _tmp8;
  const Scalar _tmp10 = -Scalar(0.010999999999999999) * _tmp9;
  const Scalar _tmp11 = 2 * _tmp0;
  const Scalar _tmp12 = _tmp11 * _tmp5;
  const Scalar _tmp13 = _tmp2 * _tmp7;
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp10 + _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp4;
  const Scalar _tmp17 = _tmp16 + position_vector(1, 0);
  const Scalar _tmp18 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp19 = Scalar(0.20999999999999999) * _tmp18 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp20 = _tmp11 * _tmp2;
  const Scalar _tmp21 = _tmp5 * _tmp7;
  const Scalar _tmp22 = _tmp20 + _tmp21;
  const Scalar _tmp23 = -Scalar(0.010999999999999999) * _tmp22;
  const Scalar _tmp24 = Scalar(0.20999999999999999) * _tmp12 - Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp25 = _tmp23 + _tmp24;
  const Scalar _tmp26 = _tmp19 + _tmp25;
  const Scalar _tmp27 = _tmp26 + position_vector(0, 0);
  const Scalar _tmp28 = Scalar(1.0) / (fh1);
  const Scalar _tmp29 = _tmp28 * fv1;
  const Scalar _tmp30 = std::asinh(_tmp29);
  const Scalar _tmp31 = Scalar(1.4083112389913199) * _tmp30;
  const Scalar _tmp32 =
      -_tmp31 * fh1 - std::sqrt(Scalar(std::pow(Scalar(-_tmp17 + p_c(1, 0)), Scalar(2)) +
                                       std::pow(Scalar(-_tmp27 + p_c(0, 0)), Scalar(2))));
  const Scalar _tmp33 = Scalar(0.71007031138673404) * _tmp28;
  const Scalar _tmp34 = _tmp32 * _tmp33;
  const Scalar _tmp35 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp36 =
      std::pow(Scalar(_tmp35 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp37 = Scalar(1.0) * _tmp30;
  const Scalar _tmp38 = _tmp27 - p_c(0, 0);
  const Scalar _tmp39 = _tmp17 - p_c(1, 0);
  const Scalar _tmp40 = std::pow(Scalar(std::pow(_tmp38, Scalar(2)) + std::pow(_tmp39, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp41 = _tmp38 * _tmp40;
  const Scalar _tmp42 = Scalar(0.20999999999999999) * _tmp20 - Scalar(0.20999999999999999) * _tmp21;
  const Scalar _tmp43 = -Scalar(0.010999999999999999) * _tmp1 -
                        Scalar(0.010999999999999999) * _tmp18 + Scalar(-0.010999999999999999);
  const Scalar _tmp44 = Scalar(0.20999999999999999) * _tmp6 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp45 = _tmp43 + _tmp44;
  const Scalar _tmp46 = _tmp42 + _tmp45;
  const Scalar _tmp47 = _tmp46 * fh1;
  const Scalar _tmp48 = Scalar(40.024799999999999) * _tmp22 + _tmp26 * fv1 + _tmp41 * _tmp47;
  const Scalar _tmp49 = -_tmp19;
  const Scalar _tmp50 = _tmp25 + _tmp49;
  const Scalar _tmp51 = -_tmp4;
  const Scalar _tmp52 = _tmp10 - _tmp14;
  const Scalar _tmp53 = _tmp51 + _tmp52;
  const Scalar _tmp54 = Scalar(1.0) * _tmp53;
  const Scalar _tmp55 = -_tmp54;
  const Scalar _tmp56 = _tmp4 + _tmp52;
  const Scalar _tmp57 = _tmp55 + _tmp56;
  const Scalar _tmp58 = _tmp15 + _tmp51;
  const Scalar _tmp59 = Scalar(1.0) / (_tmp55 + _tmp58);
  const Scalar _tmp60 = _tmp23 - _tmp24;
  const Scalar _tmp61 = _tmp19 + _tmp60;
  const Scalar _tmp62 = _tmp49 + _tmp60;
  const Scalar _tmp63 = Scalar(1.0) * _tmp62;
  const Scalar _tmp64 = _tmp59 * (-_tmp61 + _tmp63);
  const Scalar _tmp65 = _tmp57 * _tmp64;
  const Scalar _tmp66 = Scalar(1.0) / (-_tmp50 + _tmp63 - _tmp65);
  const Scalar _tmp67 = Scalar(1.0) * _tmp66;
  const Scalar _tmp68 = _tmp57 * _tmp59;
  const Scalar _tmp69 = _tmp67 * _tmp68;
  const Scalar _tmp70 = -Scalar(1.0) * _tmp67 + Scalar(1.0) * _tmp69;
  const Scalar _tmp71 = _tmp61 + position_vector(0, 0);
  const Scalar _tmp72 = _tmp71 - p_b(0, 0);
  const Scalar _tmp73 = _tmp58 + position_vector(1, 0);
  const Scalar _tmp74 = _tmp73 - p_b(1, 0);
  const Scalar _tmp75 = std::pow(Scalar(std::pow(_tmp72, Scalar(2)) + std::pow(_tmp74, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp76 = _tmp74 * _tmp75;
  const Scalar _tmp77 = _tmp72 * _tmp75;
  const Scalar _tmp78 = _tmp62 + position_vector(0, 0);
  const Scalar _tmp79 = _tmp78 - p_a(0, 0);
  const Scalar _tmp80 = _tmp53 + position_vector(1, 0);
  const Scalar _tmp81 = _tmp80 - p_a(1, 0);
  const Scalar _tmp82 =
      std::sqrt(Scalar(std::pow(_tmp79, Scalar(2)) + std::pow(_tmp81, Scalar(2))));
  const Scalar _tmp83 = Scalar(1.0) / (_tmp82);
  const Scalar _tmp84 = Scalar(1.0) / (_tmp79);
  const Scalar _tmp85 = _tmp82 * _tmp84;
  const Scalar _tmp86 = _tmp85 * (-_tmp53 * _tmp79 * _tmp83 + _tmp62 * _tmp81 * _tmp83);
  const Scalar _tmp87 = _tmp58 * _tmp77 - _tmp61 * _tmp76 + _tmp77 * _tmp86;
  const Scalar _tmp88 = _tmp81 * _tmp84;
  const Scalar _tmp89 = Scalar(1.0) / (-_tmp76 + _tmp77 * _tmp88);
  const Scalar _tmp90 = Scalar(1.0) * _tmp89;
  const Scalar _tmp91 = _tmp43 - _tmp44;
  const Scalar _tmp92 = _tmp42 + _tmp91;
  const Scalar _tmp93 = -_tmp42;
  const Scalar _tmp94 = _tmp91 + _tmp93;
  const Scalar _tmp95 = _tmp77 * _tmp94;
  const Scalar _tmp96 = _tmp76 * _tmp92 - _tmp88 * _tmp95;
  const Scalar _tmp97 = _tmp90 * _tmp96;
  const Scalar _tmp98 = -_tmp77 * _tmp92 + _tmp95;
  const Scalar _tmp99 = _tmp64 * _tmp97 - _tmp90 * _tmp98;
  const Scalar _tmp100 = _tmp50 + position_vector(0, 0);
  const Scalar _tmp101 = _tmp100 - p_d(0, 0);
  const Scalar _tmp102 = _tmp56 + position_vector(1, 0);
  const Scalar _tmp103 = _tmp102 - p_d(1, 0);
  const Scalar _tmp104 =
      std::pow(Scalar(std::pow(_tmp101, Scalar(2)) + std::pow(_tmp103, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp105 = _tmp101 * _tmp104;
  const Scalar _tmp106 = _tmp103 * _tmp104;
  const Scalar _tmp107 = _tmp105 * _tmp88 - _tmp106;
  const Scalar _tmp108 = _tmp107 * _tmp89;
  const Scalar _tmp109 = _tmp45 + _tmp93;
  const Scalar _tmp110 = _tmp88 * _tmp94;
  const Scalar _tmp111 = -_tmp105 * _tmp110 + _tmp106 * _tmp109 - _tmp108 * _tmp96;
  const Scalar _tmp112 =
      -_tmp105 * _tmp109 + _tmp105 * _tmp94 - _tmp108 * _tmp98 - _tmp111 * _tmp64;
  const Scalar _tmp113 = Scalar(1.0) / (_tmp112);
  const Scalar _tmp114 = _tmp105 * _tmp56 + _tmp105 * _tmp86 - _tmp106 * _tmp50 - _tmp108 * _tmp87;
  const Scalar _tmp115 = _tmp113 * _tmp114;
  const Scalar _tmp116 = Scalar(1.0) / (_tmp114);
  const Scalar _tmp117 = _tmp112 * _tmp116;
  const Scalar _tmp118 = _tmp117 * (-_tmp115 * _tmp99 - _tmp87 * _tmp90);
  const Scalar _tmp119 = _tmp57 * _tmp66;
  const Scalar _tmp120 = _tmp113 * (_tmp118 + _tmp99);
  const Scalar _tmp121 = -_tmp111 * _tmp120 + _tmp118 * _tmp119 - _tmp97;
  const Scalar _tmp122 = Scalar(1.0) * _tmp59;
  const Scalar _tmp123 = _tmp39 * _tmp40;
  const Scalar _tmp124 = Scalar(1.0) * _tmp123 * (_tmp118 * _tmp67 - _tmp121 * _tmp122);
  const Scalar _tmp125 = -_tmp123 * _tmp47 - _tmp16 * fv1 - Scalar(40.024799999999999) * _tmp9;
  const Scalar _tmp126 = _tmp64 * _tmp67;
  const Scalar _tmp127 = _tmp65 * _tmp67 + Scalar(1.0);
  const Scalar _tmp128 = -Scalar(1.0) * _tmp122 * _tmp127 + Scalar(1.0) * _tmp126;
  const Scalar _tmp129 = _tmp117 * _tmp67;
  const Scalar _tmp130 = Scalar(1.0) * _tmp116;
  const Scalar _tmp131 = -_tmp111 * _tmp130 + _tmp129 * _tmp57;
  const Scalar _tmp132 = _tmp123 * _tmp26 - _tmp16 * _tmp41;
  const Scalar _tmp133 = Scalar(1.0) * _tmp132;
  const Scalar _tmp134 = _tmp133 * (-_tmp122 * _tmp131 + _tmp129);
  const Scalar _tmp135 = _tmp88 * _tmp89;
  const Scalar _tmp136 = _tmp110 + _tmp135 * _tmp96;
  const Scalar _tmp137 = _tmp135 * _tmp98 - _tmp136 * _tmp64 - _tmp94;
  const Scalar _tmp138 = _tmp117 * (-_tmp115 * _tmp137 + _tmp135 * _tmp87 - _tmp86);
  const Scalar _tmp139 = _tmp113 * (_tmp137 + _tmp138);
  const Scalar _tmp140 = -_tmp111 * _tmp139 + _tmp119 * _tmp138 + _tmp136;
  const Scalar _tmp141 = Scalar(1.0) * _tmp41 * (-_tmp122 * _tmp140 + _tmp138 * _tmp67);
  const Scalar _tmp142 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp143 = _tmp54 * _tmp64 + _tmp63;
  const Scalar _tmp144 = 0;
  const Scalar _tmp145 = _tmp143 * _tmp66;
  const Scalar _tmp146 = -_tmp111 * _tmp144 - _tmp145 * _tmp57 + _tmp55;
  const Scalar _tmp147 =
      _tmp124 * fh1 + _tmp125 * _tmp128 + _tmp134 * fh1 + _tmp141 * fh1 +
      Scalar(1.0) * _tmp142 * (-_tmp122 * _tmp146 - _tmp143 * _tmp67 + Scalar(1.0)) +
      _tmp48 * _tmp70;
  const Scalar _tmp148 = _tmp89 * (-_tmp107 * _tmp139 - _tmp88);
  const Scalar _tmp149 = _tmp41 * _tmp85 * (_tmp105 * _tmp139 + _tmp148 * _tmp77 + Scalar(1.0));
  const Scalar _tmp150 =
      _tmp132 * _tmp85 * (_tmp105 * _tmp130 - _tmp107 * _tmp116 * _tmp77 * _tmp90);
  const Scalar _tmp151 = _tmp89 * (-_tmp107 * _tmp120 + Scalar(1.0));
  const Scalar _tmp152 = _tmp123 * _tmp85 * (_tmp105 * _tmp120 + _tmp151 * _tmp77);
  const Scalar _tmp153 = -_tmp142 * _tmp85 * (_tmp105 * _tmp144 - _tmp108 * _tmp144 * _tmp77) -
                         _tmp149 * fh1 - _tmp150 * fh1 - _tmp152 * fh1;
  const Scalar _tmp154 = Scalar(1.0) / (_tmp153);
  const Scalar _tmp155 = std::asinh(_tmp147 * _tmp154);
  const Scalar _tmp156 = Scalar(1.0) * _tmp155;
  const Scalar _tmp157 = Scalar(1.4083112389913199) * _tmp153;
  const Scalar _tmp158 =
      -_tmp155 * _tmp157 - std::sqrt(Scalar(std::pow(Scalar(-_tmp78 + p_a(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp80 + p_a(1, 0)), Scalar(2))));
  const Scalar _tmp159 = Scalar(0.71007031138673404) * _tmp154;
  const Scalar _tmp160 = _tmp158 * _tmp159;
  const Scalar _tmp161 = -_tmp149 - _tmp150 - _tmp152;
  const Scalar _tmp162 = Scalar(1.4083112389913199) * _tmp161;
  const Scalar _tmp163 = _tmp41 * _tmp46;
  const Scalar _tmp164 = _tmp123 * _tmp46;
  const Scalar _tmp165 = std::pow(_tmp153, Scalar(-2));
  const Scalar _tmp166 = _tmp161 * _tmp165;
  const Scalar _tmp167 = (-_tmp147 * _tmp166 + _tmp154 * (_tmp124 - _tmp128 * _tmp164 + _tmp134 +
                                                          _tmp141 + _tmp163 * _tmp70)) /
                         std::sqrt(Scalar(std::pow(_tmp147, Scalar(2)) * _tmp165 + 1));
  const Scalar _tmp168 = _tmp127 * _tmp59;
  const Scalar _tmp169 = _tmp140 * _tmp41 * _tmp59;
  const Scalar _tmp170 = _tmp131 * _tmp132 * _tmp59;
  const Scalar _tmp171 = _tmp121 * _tmp123 * _tmp59;
  const Scalar _tmp172 = _tmp48 * _tmp67;
  const Scalar _tmp173 = _tmp125 * _tmp168 + _tmp142 * _tmp146 * _tmp59 + _tmp169 * fh1 +
                         _tmp170 * fh1 + _tmp171 * fh1 - _tmp172 * _tmp68;
  const Scalar _tmp174 = _tmp123 * _tmp151;
  const Scalar _tmp175 = _tmp142 * _tmp144;
  const Scalar _tmp176 = _tmp116 * _tmp133;
  const Scalar _tmp177 = _tmp176 * fh1;
  const Scalar _tmp178 = _tmp148 * _tmp41;
  const Scalar _tmp179 = -_tmp108 * _tmp175 - _tmp108 * _tmp177 + _tmp174 * fh1 + _tmp178 * fh1;
  const Scalar _tmp180 = Scalar(1.0) / (_tmp179);
  const Scalar _tmp181 = std::asinh(_tmp173 * _tmp180);
  const Scalar _tmp182 = Scalar(1.4083112389913199) * _tmp179;
  const Scalar _tmp183 =
      -_tmp181 * _tmp182 - std::sqrt(Scalar(std::pow(Scalar(-_tmp71 + p_b(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp73 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp184 = Scalar(0.71007031138673404) * _tmp180;
  const Scalar _tmp185 = _tmp183 * _tmp184;
  const Scalar _tmp186 = Scalar(1.0) * _tmp181;
  const Scalar _tmp187 = -_tmp108 * _tmp176 + _tmp174 + _tmp178;
  const Scalar _tmp188 = Scalar(1.4083112389913199) * _tmp187;
  const Scalar _tmp189 = std::pow(_tmp179, Scalar(-2));
  const Scalar _tmp190 = _tmp187 * _tmp189;
  const Scalar _tmp191 = (-_tmp173 * _tmp190 + _tmp180 * (-_tmp163 * _tmp69 - _tmp164 * _tmp168 +
                                                          _tmp169 + _tmp170 + _tmp171)) /
                         std::sqrt(Scalar(std::pow(_tmp173, Scalar(2)) * _tmp189 + 1));
  const Scalar _tmp192 = _tmp138 * _tmp41 * _tmp66;
  const Scalar _tmp193 = _tmp129 * _tmp132;
  const Scalar _tmp194 = _tmp118 * _tmp123 * _tmp66;
  const Scalar _tmp195 = -_tmp125 * _tmp126 + _tmp142 * _tmp145 + _tmp172 - _tmp192 * fh1 -
                         _tmp193 * fh1 - _tmp194 * fh1;
  const Scalar _tmp196 = _tmp120 * _tmp123;
  const Scalar _tmp197 = _tmp139 * _tmp41;
  const Scalar _tmp198 = _tmp175 + _tmp177 + _tmp196 * fh1 + _tmp197 * fh1;
  const Scalar _tmp199 = Scalar(1.0) / (_tmp198);
  const Scalar _tmp200 = std::asinh(_tmp195 * _tmp199);
  const Scalar _tmp201 = Scalar(1.0) * _tmp200;
  const Scalar _tmp202 = Scalar(1.4083112389913199) * _tmp198;
  const Scalar _tmp203 =
      -_tmp200 * _tmp202 - std::sqrt(Scalar(std::pow(Scalar(-_tmp100 + p_d(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp102 + p_d(1, 0)), Scalar(2))));
  const Scalar _tmp204 = Scalar(0.71007031138673404) * _tmp199;
  const Scalar _tmp205 = _tmp203 * _tmp204;
  const Scalar _tmp206 = _tmp176 + _tmp196 + _tmp197;
  const Scalar _tmp207 = Scalar(1.4083112389913199) * _tmp206;
  const Scalar _tmp208 = std::pow(_tmp198, Scalar(-2));
  const Scalar _tmp209 = _tmp206 * _tmp208;
  const Scalar _tmp210 = (-_tmp195 * _tmp209 + _tmp199 * (_tmp126 * _tmp164 + _tmp163 * _tmp67 -
                                                          _tmp192 - _tmp193 - _tmp194)) /
                         std::sqrt(Scalar(std::pow(_tmp195, Scalar(2)) * _tmp208 + 1));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = Scalar(1.4083112389913199) * fh1 *
                   (Scalar(1.0) * _tmp35 * _tmp36 * fv1 * std::cosh(_tmp37) -
                    (-Scalar(0.71007031138673404) * _tmp32 * _tmp35 +
                     _tmp33 * (Scalar(1.4083112389913199) * _tmp29 * _tmp36 - _tmp31)) *
                        std::cosh(_tmp34)) -
               Scalar(1.4083112389913199) * std::sinh(_tmp34) -
               Scalar(1.4083112389913199) * std::sinh(_tmp37);
  _res(1, 0) = _tmp157 * (-Scalar(1.0) * _tmp167 * std::cosh(_tmp156) -
                          (-Scalar(0.71007031138673404) * _tmp158 * _tmp166 +
                           _tmp159 * (-_tmp155 * _tmp162 - _tmp157 * _tmp167)) *
                              std::cosh(_tmp160)) +
               _tmp162 * (-std::sinh(_tmp156) - std::sinh(_tmp160));
  _res(2, 0) = _tmp182 * (-Scalar(1.0) * _tmp191 * std::cosh(_tmp186) -
                          (-Scalar(0.71007031138673404) * _tmp183 * _tmp190 +
                           _tmp184 * (-_tmp181 * _tmp188 - _tmp182 * _tmp191)) *
                              std::cosh(_tmp185)) +
               _tmp188 * (-std::sinh(_tmp185) - std::sinh(_tmp186));
  _res(3, 0) = _tmp202 * (-Scalar(1.0) * _tmp210 * std::cosh(_tmp201) -
                          (-Scalar(0.71007031138673404) * _tmp203 * _tmp209 +
                           _tmp204 * (-_tmp200 * _tmp207 - _tmp202 * _tmp210)) *
                              std::cosh(_tmp205)) +
               _tmp207 * (-std::sinh(_tmp201) - std::sinh(_tmp205));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
