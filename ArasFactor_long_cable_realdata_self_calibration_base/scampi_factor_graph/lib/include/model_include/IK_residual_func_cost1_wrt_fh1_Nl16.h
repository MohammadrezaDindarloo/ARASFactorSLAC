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
 * Symbolic function: IK_residual_func_cost1_wrt_fh1_Nl16
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1WrtFh1Nl16(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 635

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (216)
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
  const Scalar _tmp8 = _tmp2 * _tmp4;
  const Scalar _tmp9 = _tmp0 * _tmp5;
  const Scalar _tmp10 = _tmp8 - _tmp9;
  const Scalar _tmp11 = -Scalar(0.010999999999999999) * _tmp10;
  const Scalar _tmp12 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp13 = 1 - 2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp12 + Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = _tmp11 + _tmp14;
  const Scalar _tmp16 = _tmp15 + _tmp7;
  const Scalar _tmp17 = _tmp16 + position_vector(1, 0);
  const Scalar _tmp18 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp19 = Scalar(0.20999999999999999) * _tmp12 +
                        Scalar(0.20999999999999999) * _tmp18 + Scalar(0.20999999999999999);
  const Scalar _tmp20 = 2 * _tmp0 * _tmp4;
  const Scalar _tmp21 = _tmp1 * _tmp5;
  const Scalar _tmp22 = _tmp20 + _tmp21;
  const Scalar _tmp23 = -Scalar(0.010999999999999999) * _tmp22;
  const Scalar _tmp24 = Scalar(0.20999999999999999) * _tmp3 - Scalar(0.20999999999999999) * _tmp6;
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
  const Scalar _tmp36 = Scalar(0.71007031138673404) * _tmp35;
  const Scalar _tmp37 =
      std::pow(Scalar(_tmp35 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp38 = Scalar(1.0) * _tmp30;
  const Scalar _tmp39 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp40 = -_tmp19;
  const Scalar _tmp41 = _tmp25 + _tmp40;
  const Scalar _tmp42 = Scalar(1.0) * _tmp41;
  const Scalar _tmp43 = -_tmp7;
  const Scalar _tmp44 = _tmp15 + _tmp43;
  const Scalar _tmp45 = Scalar(1.0) * _tmp44;
  const Scalar _tmp46 = -_tmp45;
  const Scalar _tmp47 = _tmp11 - _tmp14;
  const Scalar _tmp48 = _tmp43 + _tmp47;
  const Scalar _tmp49 = Scalar(1.0) / (_tmp46 + _tmp48);
  const Scalar _tmp50 = _tmp23 - _tmp24;
  const Scalar _tmp51 = _tmp40 + _tmp50;
  const Scalar _tmp52 = _tmp49 * (_tmp42 - _tmp51);
  const Scalar _tmp53 = _tmp42 + _tmp45 * _tmp52;
  const Scalar _tmp54 = Scalar(0.20999999999999999) * _tmp20 - Scalar(0.20999999999999999) * _tmp21;
  const Scalar _tmp55 = -_tmp54;
  const Scalar _tmp56 =
      -Scalar(0.010999999999999999) * _tmp13 - Scalar(0.010999999999999999) * _tmp18;
  const Scalar _tmp57 = Scalar(0.20999999999999999) * _tmp8 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp58 = _tmp56 + _tmp57;
  const Scalar _tmp59 = _tmp55 + _tmp58;
  const Scalar _tmp60 = _tmp19 + _tmp50;
  const Scalar _tmp61 = _tmp60 + position_vector(0, 0);
  const Scalar _tmp62 = _tmp61 - p_b(0, 0);
  const Scalar _tmp63 = _tmp47 + _tmp7;
  const Scalar _tmp64 = _tmp63 + position_vector(1, 0);
  const Scalar _tmp65 = _tmp64 - p_b(1, 0);
  const Scalar _tmp66 = std::pow(Scalar(std::pow(_tmp62, Scalar(2)) + std::pow(_tmp65, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp67 = _tmp62 * _tmp66;
  const Scalar _tmp68 = _tmp59 * _tmp67;
  const Scalar _tmp69 = _tmp56 - _tmp57;
  const Scalar _tmp70 = _tmp54 + _tmp69;
  const Scalar _tmp71 = _tmp65 * _tmp66;
  const Scalar _tmp72 = _tmp41 + position_vector(0, 0);
  const Scalar _tmp73 = _tmp72 - p_d(0, 0);
  const Scalar _tmp74 = Scalar(1.0) / (_tmp73);
  const Scalar _tmp75 = _tmp44 + position_vector(1, 0);
  const Scalar _tmp76 = _tmp75 - p_d(1, 0);
  const Scalar _tmp77 = _tmp74 * _tmp76;
  const Scalar _tmp78 = _tmp51 + position_vector(0, 0);
  const Scalar _tmp79 = _tmp78 - p_a(0, 0);
  const Scalar _tmp80 = _tmp48 + position_vector(1, 0);
  const Scalar _tmp81 = _tmp80 - p_a(1, 0);
  const Scalar _tmp82 = std::pow(Scalar(std::pow(_tmp79, Scalar(2)) + std::pow(_tmp81, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp83 = _tmp79 * _tmp82;
  const Scalar _tmp84 = _tmp59 * _tmp83;
  const Scalar _tmp85 = _tmp55 + _tmp69;
  const Scalar _tmp86 = _tmp81 * _tmp82;
  const Scalar _tmp87 = -_tmp77 * _tmp84 + _tmp85 * _tmp86;
  const Scalar _tmp88 = Scalar(1.0) / (_tmp77 * _tmp83 - _tmp86);
  const Scalar _tmp89 = _tmp67 * _tmp77 - _tmp71;
  const Scalar _tmp90 = _tmp88 * _tmp89;
  const Scalar _tmp91 = -_tmp68 * _tmp77 + _tmp70 * _tmp71 - _tmp87 * _tmp90;
  const Scalar _tmp92 = -_tmp83 * _tmp85 + _tmp84;
  const Scalar _tmp93 = -_tmp52 * _tmp91 - _tmp67 * _tmp70 + _tmp68 - _tmp90 * _tmp92;
  const Scalar _tmp94 = Scalar(1.0) / (_tmp93);
  const Scalar _tmp95 = 0;
  const Scalar _tmp96 =
      std::sqrt(Scalar(std::pow(_tmp73, Scalar(2)) + std::pow(_tmp76, Scalar(2))));
  const Scalar _tmp97 = _tmp74 * _tmp96;
  const Scalar _tmp98 = _tmp27 - p_c(0, 0);
  const Scalar _tmp99 = _tmp17 - p_c(1, 0);
  const Scalar _tmp100 = std::pow(Scalar(std::pow(_tmp98, Scalar(2)) + std::pow(_tmp99, Scalar(2))),
                                  Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp101 = _tmp100 * _tmp98;
  const Scalar _tmp102 = _tmp100 * _tmp99;
  const Scalar _tmp103 = -_tmp101 * _tmp16 + _tmp102 * _tmp26;
  const Scalar _tmp104 = Scalar(1.0) / (_tmp96);
  const Scalar _tmp105 = _tmp97 * (_tmp104 * _tmp41 * _tmp76 - _tmp104 * _tmp44 * _tmp73);
  const Scalar _tmp106 = _tmp105 * _tmp83 + _tmp48 * _tmp83 - _tmp51 * _tmp86;
  const Scalar _tmp107 = _tmp105 * _tmp67 - _tmp106 * _tmp90 - _tmp60 * _tmp71 + _tmp63 * _tmp67;
  const Scalar _tmp108 = Scalar(1.0) / (_tmp107);
  const Scalar _tmp109 = Scalar(1.0) * _tmp108;
  const Scalar _tmp110 = Scalar(1.0) * _tmp88;
  const Scalar _tmp111 = _tmp108 * _tmp110 * _tmp89;
  const Scalar _tmp112 = _tmp103 * _tmp97 * (_tmp109 * _tmp67 - _tmp111 * _tmp83);
  const Scalar _tmp113 = _tmp77 * _tmp88;
  const Scalar _tmp114 = _tmp113 * _tmp87 + _tmp59 * _tmp77;
  const Scalar _tmp115 = _tmp113 * _tmp92 - _tmp114 * _tmp52 - _tmp59;
  const Scalar _tmp116 = _tmp107 * _tmp94;
  const Scalar _tmp117 = _tmp108 * _tmp93;
  const Scalar _tmp118 = _tmp117 * (-_tmp105 + _tmp106 * _tmp113 - _tmp115 * _tmp116);
  const Scalar _tmp119 = _tmp94 * (_tmp115 + _tmp118);
  const Scalar _tmp120 = -_tmp119 * _tmp89 - _tmp77;
  const Scalar _tmp121 = _tmp83 * _tmp88;
  const Scalar _tmp122 = _tmp101 * _tmp97 * (_tmp119 * _tmp67 + _tmp120 * _tmp121 + Scalar(1.0));
  const Scalar _tmp123 = _tmp110 * _tmp87;
  const Scalar _tmp124 = -_tmp110 * _tmp92 + _tmp123 * _tmp52;
  const Scalar _tmp125 = _tmp117 * (-_tmp106 * _tmp110 - _tmp116 * _tmp124);
  const Scalar _tmp126 = _tmp94 * (_tmp124 + _tmp125);
  const Scalar _tmp127 = -_tmp126 * _tmp89 + Scalar(1.0);
  const Scalar _tmp128 = _tmp102 * _tmp97 * (_tmp121 * _tmp127 + _tmp126 * _tmp67);
  const Scalar _tmp129 = -_tmp112 * fh1 - _tmp122 * fh1 - _tmp128 * fh1 -
                         _tmp39 * _tmp97 * (_tmp67 * _tmp95 - _tmp83 * _tmp90 * _tmp95);
  const Scalar _tmp130 = std::pow(_tmp129, Scalar(-2));
  const Scalar _tmp131 = -_tmp112 - _tmp122 - _tmp128;
  const Scalar _tmp132 = _tmp130 * _tmp131;
  const Scalar _tmp133 = Scalar(0.71007031138673404) * _tmp132;
  const Scalar _tmp134 = Scalar(1.0) / (_tmp129);
  const Scalar _tmp135 = _tmp46 + _tmp63;
  const Scalar _tmp136 = _tmp135 * _tmp52;
  const Scalar _tmp137 = Scalar(1.0) / (-_tmp136 + _tmp42 - _tmp60);
  const Scalar _tmp138 = Scalar(1.0) * _tmp137;
  const Scalar _tmp139 = _tmp135 * _tmp137;
  const Scalar _tmp140 = _tmp114 + _tmp118 * _tmp139 - _tmp119 * _tmp91;
  const Scalar _tmp141 = Scalar(1.0) * _tmp49;
  const Scalar _tmp142 = Scalar(1.0) * _tmp101 * (_tmp118 * _tmp138 - _tmp140 * _tmp141);
  const Scalar _tmp143 = _tmp135 * _tmp138;
  const Scalar _tmp144 = -_tmp109 * _tmp91 + _tmp117 * _tmp143;
  const Scalar _tmp145 = _tmp117 * _tmp138;
  const Scalar _tmp146 = Scalar(1.0) * _tmp103;
  const Scalar _tmp147 = _tmp146 * (-_tmp141 * _tmp144 + _tmp145);
  const Scalar _tmp148 = -_tmp123 + _tmp125 * _tmp139 - _tmp126 * _tmp91;
  const Scalar _tmp149 = Scalar(1.0) * _tmp102 * (_tmp125 * _tmp138 - _tmp141 * _tmp148);
  const Scalar _tmp150 = _tmp54 + _tmp58;
  const Scalar _tmp151 = _tmp150 * fh1;
  const Scalar _tmp152 = _tmp101 * _tmp151 + Scalar(40.024799999999999) * _tmp22 + _tmp26 * fv1;
  const Scalar _tmp153 = _tmp143 * _tmp49;
  const Scalar _tmp154 = -Scalar(1.0) * _tmp138 + Scalar(1.0) * _tmp153;
  const Scalar _tmp155 = -Scalar(40.024799999999999) * _tmp10 - _tmp102 * _tmp151 - _tmp16 * fv1;
  const Scalar _tmp156 = _tmp49 * (_tmp136 * _tmp138 + Scalar(1.0));
  const Scalar _tmp157 = _tmp138 * _tmp52;
  const Scalar _tmp158 = -Scalar(1.0) * _tmp156 + Scalar(1.0) * _tmp157;
  const Scalar _tmp159 = _tmp137 * _tmp53;
  const Scalar _tmp160 = -_tmp135 * _tmp159 + _tmp46 - _tmp91 * _tmp95;
  const Scalar _tmp161 =
      _tmp142 * fh1 + _tmp147 * fh1 + _tmp149 * fh1 + _tmp152 * _tmp154 + _tmp155 * _tmp158 +
      Scalar(1.0) * _tmp39 * (-_tmp138 * _tmp53 - _tmp141 * _tmp160 + Scalar(1.0));
  const Scalar _tmp162 = std::asinh(_tmp134 * _tmp161);
  const Scalar _tmp163 = Scalar(1.4083112389913199) * _tmp162;
  const Scalar _tmp164 = Scalar(1.4083112389913199) * _tmp129;
  const Scalar _tmp165 = _tmp102 * _tmp150;
  const Scalar _tmp166 = _tmp101 * _tmp150;
  const Scalar _tmp167 = (-_tmp132 * _tmp161 + _tmp134 * (_tmp142 + _tmp147 + _tmp149 +
                                                          _tmp154 * _tmp166 - _tmp158 * _tmp165)) /
                         std::sqrt(Scalar(_tmp130 * std::pow(_tmp161, Scalar(2)) + 1));
  const Scalar _tmp168 = Scalar(0.71007031138673404) * _tmp134;
  const Scalar _tmp169 =
      -_tmp129 * _tmp163 - std::sqrt(Scalar(std::pow(Scalar(-_tmp72 + p_d(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp75 + p_d(1, 0)), Scalar(2))));
  const Scalar _tmp170 = _tmp168 * _tmp169;
  const Scalar _tmp171 = Scalar(1.0) * _tmp162;
  const Scalar _tmp172 = _tmp103 * _tmp144 * _tmp49;
  const Scalar _tmp173 = _tmp102 * _tmp148 * _tmp49;
  const Scalar _tmp174 = _tmp101 * _tmp140 * _tmp49;
  const Scalar _tmp175 = _tmp138 * _tmp152;
  const Scalar _tmp176 = -_tmp135 * _tmp175 * _tmp49 + _tmp155 * _tmp156 +
                         _tmp160 * _tmp39 * _tmp49 + _tmp172 * fh1 + _tmp173 * fh1 + _tmp174 * fh1;
  const Scalar _tmp177 = _tmp102 * _tmp127 * _tmp88;
  const Scalar _tmp178 = _tmp101 * _tmp120 * _tmp88;
  const Scalar _tmp179 = _tmp39 * _tmp95;
  const Scalar _tmp180 = _tmp103 * _tmp111;
  const Scalar _tmp181 = _tmp177 * fh1 + _tmp178 * fh1 - _tmp179 * _tmp90 - _tmp180 * fh1;
  const Scalar _tmp182 = Scalar(1.0) / (_tmp181);
  const Scalar _tmp183 = std::asinh(_tmp176 * _tmp182);
  const Scalar _tmp184 = Scalar(1.4083112389913199) * _tmp181;
  const Scalar _tmp185 =
      -_tmp183 * _tmp184 - std::sqrt(Scalar(std::pow(Scalar(-_tmp78 + p_a(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp80 + p_a(1, 0)), Scalar(2))));
  const Scalar _tmp186 = Scalar(0.71007031138673404) * _tmp182;
  const Scalar _tmp187 = _tmp185 * _tmp186;
  const Scalar _tmp188 = Scalar(1.0) * _tmp183;
  const Scalar _tmp189 = _tmp177 + _tmp178 - _tmp180;
  const Scalar _tmp190 = Scalar(1.4083112389913199) * _tmp189;
  const Scalar _tmp191 = std::pow(_tmp181, Scalar(-2));
  const Scalar _tmp192 = _tmp189 * _tmp191;
  const Scalar _tmp193 = Scalar(0.71007031138673404) * _tmp192;
  const Scalar _tmp194 = (-_tmp176 * _tmp192 + _tmp182 * (-_tmp153 * _tmp166 - _tmp156 * _tmp165 +
                                                          _tmp172 + _tmp173 + _tmp174)) /
                         std::sqrt(Scalar(std::pow(_tmp176, Scalar(2)) * _tmp191 + 1));
  const Scalar _tmp195 = _tmp101 * _tmp119;
  const Scalar _tmp196 = _tmp108 * _tmp146;
  const Scalar _tmp197 = _tmp102 * _tmp126;
  const Scalar _tmp198 = _tmp179 + _tmp195 * fh1 + _tmp196 * fh1 + _tmp197 * fh1;
  const Scalar _tmp199 = Scalar(1.0) / (_tmp198);
  const Scalar _tmp200 = _tmp103 * _tmp145;
  const Scalar _tmp201 = _tmp102 * _tmp125 * _tmp137;
  const Scalar _tmp202 = _tmp101 * _tmp118 * _tmp137;
  const Scalar _tmp203 = -_tmp155 * _tmp157 + _tmp159 * _tmp39 + _tmp175 - _tmp200 * fh1 -
                         _tmp201 * fh1 - _tmp202 * fh1;
  const Scalar _tmp204 = std::asinh(_tmp199 * _tmp203);
  const Scalar _tmp205 = Scalar(1.4083112389913199) * _tmp198;
  const Scalar _tmp206 =
      -_tmp204 * _tmp205 - std::sqrt(Scalar(std::pow(Scalar(-_tmp61 + p_b(0, 0)), Scalar(2)) +
                                            std::pow(Scalar(-_tmp64 + p_b(1, 0)), Scalar(2))));
  const Scalar _tmp207 = Scalar(0.71007031138673404) * _tmp199;
  const Scalar _tmp208 = _tmp206 * _tmp207;
  const Scalar _tmp209 = Scalar(1.0) * _tmp204;
  const Scalar _tmp210 = _tmp195 + _tmp196 + _tmp197;
  const Scalar _tmp211 = Scalar(1.4083112389913199) * _tmp210;
  const Scalar _tmp212 = std::pow(_tmp198, Scalar(-2));
  const Scalar _tmp213 = _tmp210 * _tmp212;
  const Scalar _tmp214 = Scalar(0.71007031138673404) * _tmp213;
  const Scalar _tmp215 =
      (_tmp199 * (_tmp138 * _tmp166 + _tmp157 * _tmp165 - _tmp200 - _tmp201 - _tmp202) -
       _tmp203 * _tmp213) /
      std::sqrt(Scalar(std::pow(_tmp203, Scalar(2)) * _tmp212 + 1));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      -Scalar(1.4083112389913199) * _tmp33 * p_c(2, 0) -
      Scalar(1.4083112389913199) * fh1 *
          (-Scalar(1.0) * _tmp35 * _tmp37 * fv1 * std::sinh(_tmp38) - _tmp36 * p_c(2, 0) -
           (-_tmp32 * _tmp36 + _tmp33 * (Scalar(1.4083112389913199) * _tmp29 * _tmp37 - _tmp31)) *
               std::sinh(_tmp34)) +
      Scalar(1.4083112389913199) * std::cosh(_tmp34) -
      Scalar(1.4083112389913199) * std::cosh(_tmp38);
  _res(1, 0) =
      -Scalar(1.4083112389913199) * _tmp131 *
          (_tmp168 * p_d(2, 0) - std::cosh(_tmp170) + std::cosh(_tmp171)) -
      _tmp164 * (-_tmp133 * p_d(2, 0) + Scalar(1.0) * _tmp167 * std::sinh(_tmp171) -
                 (-_tmp133 * _tmp169 + _tmp168 * (-_tmp131 * _tmp163 - _tmp164 * _tmp167)) *
                     std::sinh(_tmp170));
  _res(2, 0) =
      -_tmp184 * (-_tmp193 * p_a(2, 0) + Scalar(1.0) * _tmp194 * std::sinh(_tmp188) -
                  (-_tmp185 * _tmp193 + _tmp186 * (-_tmp183 * _tmp190 - _tmp184 * _tmp194)) *
                      std::sinh(_tmp187)) -
      _tmp190 * (_tmp186 * p_a(2, 0) - std::cosh(_tmp187) + std::cosh(_tmp188));
  _res(3, 0) =
      -_tmp205 * (-_tmp214 * p_b(2, 0) + Scalar(1.0) * _tmp215 * std::sinh(_tmp209) -
                  (-_tmp206 * _tmp214 + _tmp207 * (-_tmp204 * _tmp211 - _tmp205 * _tmp215)) *
                      std::sinh(_tmp208)) -
      _tmp211 * (_tmp207 * p_b(2, 0) - std::cosh(_tmp208) + std::cosh(_tmp209));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
