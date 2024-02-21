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
 * Symbolic function: IK_residual_func_cost2_wrt_fh1_Nl8
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2WrtFh1Nl8(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 642

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (212)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp1 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp2 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp3 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp5 = -_tmp4;
  const Scalar _tmp6 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp7 = 2 * _tmp6;
  const Scalar _tmp8 = _tmp2 * _tmp7;
  const Scalar _tmp9 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp10 = _tmp0 * _tmp9;
  const Scalar _tmp11 = -_tmp10 + _tmp8;
  const Scalar _tmp12 = -Scalar(0.010999999999999999) * _tmp11;
  const Scalar _tmp13 = _tmp0 * _tmp7;
  const Scalar _tmp14 = _tmp2 * _tmp9;
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp12 + _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp5;
  const Scalar _tmp18 = _tmp17 + position_vector(1, 0);
  const Scalar _tmp19 = -2 * std::pow(_tmp6, Scalar(2));
  const Scalar _tmp20 = Scalar(0.20999999999999999) * _tmp19 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp21 = 2 * _tmp0 * _tmp2;
  const Scalar _tmp22 = _tmp6 * _tmp9;
  const Scalar _tmp23 = _tmp21 + _tmp22;
  const Scalar _tmp24 = -Scalar(0.010999999999999999) * _tmp23;
  const Scalar _tmp25 = Scalar(0.20999999999999999) * _tmp13 - Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp26 = _tmp24 - _tmp25;
  const Scalar _tmp27 = _tmp20 + _tmp26;
  const Scalar _tmp28 = _tmp27 + position_vector(0, 0);
  const Scalar _tmp29 = Scalar(1.0) / (fh1);
  const Scalar _tmp30 = _tmp29 * fv1;
  const Scalar _tmp31 = std::asinh(_tmp30);
  const Scalar _tmp32 = Scalar(1.4083112389913199) * _tmp31;
  const Scalar _tmp33 =
      -_tmp32 * fh1 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp28), Scalar(2)) +
                     Scalar(0.77439999999999998) *
                         std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp18 - 1), Scalar(2))));
  const Scalar _tmp34 = Scalar(0.71007031138673404) * _tmp29;
  const Scalar _tmp35 = _tmp33 * _tmp34;
  const Scalar _tmp36 = Scalar(1.0) * _tmp31;
  const Scalar _tmp37 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp38 =
      std::pow(Scalar(_tmp37 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp39 = _tmp24 + _tmp25;
  const Scalar _tmp40 = _tmp20 + _tmp39;
  const Scalar _tmp41 = _tmp40 + position_vector(0, 0);
  const Scalar _tmp42 = _tmp16 + _tmp4;
  const Scalar _tmp43 = _tmp42 + position_vector(1, 0);
  const Scalar _tmp44 = -_tmp20;
  const Scalar _tmp45 = _tmp39 + _tmp44;
  const Scalar _tmp46 = _tmp45 + position_vector(0, 0);
  const Scalar _tmp47 = _tmp46 + Scalar(125.0);
  const Scalar _tmp48 = _tmp12 - _tmp15;
  const Scalar _tmp49 = _tmp4 + _tmp48;
  const Scalar _tmp50 = _tmp49 + position_vector(1, 0);
  const Scalar _tmp51 = _tmp50 + Scalar(-110.0);
  const Scalar _tmp52 = std::pow(Scalar(std::pow(_tmp47, Scalar(2)) + std::pow(_tmp51, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp53 = _tmp47 * _tmp52;
  const Scalar _tmp54 = Scalar(0.20999999999999999) * _tmp21 - Scalar(0.20999999999999999) * _tmp22;
  const Scalar _tmp55 = -_tmp54;
  const Scalar _tmp56 = -Scalar(0.010999999999999999) * _tmp1 -
                        Scalar(0.010999999999999999) * _tmp19 + Scalar(-0.010999999999999999);
  const Scalar _tmp57 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp58 = _tmp56 + _tmp57;
  const Scalar _tmp59 = _tmp55 + _tmp58;
  const Scalar _tmp60 = _tmp54 + _tmp58;
  const Scalar _tmp61 = _tmp43 + Scalar(-110.0);
  const Scalar _tmp62 = _tmp41 + Scalar(-125.0);
  const Scalar _tmp63 = Scalar(1.0) / (_tmp62);
  const Scalar _tmp64 = _tmp61 * _tmp63;
  const Scalar _tmp65 = _tmp60 * _tmp64;
  const Scalar _tmp66 = _tmp51 * _tmp52;
  const Scalar _tmp67 = _tmp53 * _tmp64 - _tmp66;
  const Scalar _tmp68 = _tmp48 + _tmp5;
  const Scalar _tmp69 = _tmp68 + position_vector(1, 0);
  const Scalar _tmp70 = _tmp69 + Scalar(110.0);
  const Scalar _tmp71 = _tmp26 + _tmp44;
  const Scalar _tmp72 = _tmp71 + position_vector(0, 0);
  const Scalar _tmp73 = _tmp72 + Scalar(125.0);
  const Scalar _tmp74 = std::pow(Scalar(std::pow(_tmp70, Scalar(2)) + std::pow(_tmp73, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp75 = _tmp73 * _tmp74;
  const Scalar _tmp76 = _tmp70 * _tmp74;
  const Scalar _tmp77 = Scalar(1.0) / (_tmp64 * _tmp75 - _tmp76);
  const Scalar _tmp78 = _tmp56 - _tmp57;
  const Scalar _tmp79 = _tmp55 + _tmp78;
  const Scalar _tmp80 = _tmp77 * (-_tmp65 * _tmp75 + _tmp76 * _tmp79);
  const Scalar _tmp81 = -_tmp53 * _tmp65 + _tmp59 * _tmp66 - _tmp67 * _tmp80;
  const Scalar _tmp82 = Scalar(1.0) * _tmp42;
  const Scalar _tmp83 = -_tmp82;
  const Scalar _tmp84 = Scalar(1.0) / (_tmp68 + _tmp83);
  const Scalar _tmp85 = Scalar(1.0) * _tmp40;
  const Scalar _tmp86 = _tmp84 * (-_tmp71 + _tmp85);
  const Scalar _tmp87 = _tmp77 * (_tmp60 * _tmp75 - _tmp75 * _tmp79);
  const Scalar _tmp88 = -_tmp53 * _tmp59 + _tmp53 * _tmp60 - _tmp67 * _tmp87 - _tmp81 * _tmp86;
  const Scalar _tmp89 = Scalar(1.0) / (_tmp88);
  const Scalar _tmp90 =
      std::sqrt(Scalar(std::pow(_tmp61, Scalar(2)) + std::pow(_tmp62, Scalar(2))));
  const Scalar _tmp91 = Scalar(1.0) / (_tmp90);
  const Scalar _tmp92 = _tmp63 * _tmp90;
  const Scalar _tmp93 = _tmp92 * (_tmp40 * _tmp61 * _tmp91 - _tmp42 * _tmp62 * _tmp91);
  const Scalar _tmp94 = _tmp77 * (_tmp68 * _tmp75 - _tmp71 * _tmp76 + _tmp75 * _tmp93);
  const Scalar _tmp95 = Scalar(1.0) * _tmp80;
  const Scalar _tmp96 = _tmp86 * _tmp95 - Scalar(1.0) * _tmp87;
  const Scalar _tmp97 = -_tmp45 * _tmp66 + _tmp49 * _tmp53 + _tmp53 * _tmp93 - _tmp67 * _tmp94;
  const Scalar _tmp98 = _tmp89 * _tmp97;
  const Scalar _tmp99 = Scalar(1.0) / (_tmp97);
  const Scalar _tmp100 = _tmp88 * _tmp99;
  const Scalar _tmp101 = _tmp100 * (-Scalar(1.0) * _tmp94 - _tmp96 * _tmp98);
  const Scalar _tmp102 = _tmp101 + _tmp96;
  const Scalar _tmp103 = _tmp102 * _tmp89;
  const Scalar _tmp104 = _tmp67 * _tmp89;
  const Scalar _tmp105 = _tmp77 * (-_tmp102 * _tmp104 + Scalar(1.0));
  const Scalar _tmp106 = _tmp18 + Scalar(110.0);
  const Scalar _tmp107 = _tmp28 + Scalar(-125.0);
  const Scalar _tmp108 =
      std::pow(Scalar(std::pow(_tmp106, Scalar(2)) + std::pow(_tmp107, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp109 = _tmp106 * _tmp108;
  const Scalar _tmp110 = _tmp109 * _tmp92 * (_tmp103 * _tmp53 + _tmp105 * _tmp75);
  const Scalar _tmp111 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp112 = _tmp82 * _tmp86 + _tmp85;
  const Scalar _tmp113 = 0;
  const Scalar _tmp114 = _tmp113 * _tmp89;
  const Scalar _tmp115 = _tmp104 * _tmp113 * _tmp77;
  const Scalar _tmp116 = _tmp107 * _tmp108;
  const Scalar _tmp117 = _tmp109 * _tmp27 - _tmp116 * _tmp17;
  const Scalar _tmp118 = Scalar(1.0) * _tmp99;
  const Scalar _tmp119 = _tmp67 * _tmp77;
  const Scalar _tmp120 = _tmp117 * _tmp92 * (-_tmp118 * _tmp119 * _tmp75 + _tmp118 * _tmp53);
  const Scalar _tmp121 = _tmp64 * _tmp80 + _tmp65;
  const Scalar _tmp122 = -_tmp121 * _tmp86 - _tmp60 + _tmp64 * _tmp87;
  const Scalar _tmp123 = _tmp100 * (-_tmp122 * _tmp98 + _tmp64 * _tmp94 - _tmp93);
  const Scalar _tmp124 = _tmp122 + _tmp123;
  const Scalar _tmp125 = _tmp124 * _tmp89;
  const Scalar _tmp126 = _tmp77 * (-_tmp104 * _tmp124 - _tmp64);
  const Scalar _tmp127 = _tmp116 * _tmp92 * (_tmp125 * _tmp53 + _tmp126 * _tmp75 + Scalar(1.0));
  const Scalar _tmp128 = -_tmp110 * fh1 - _tmp111 * _tmp92 * (_tmp114 * _tmp53 - _tmp115 * _tmp75) -
                         _tmp120 * fh1 - _tmp127 * fh1;
  const Scalar _tmp129 = Scalar(1.0) / (_tmp128);
  const Scalar _tmp130 = _tmp81 * _tmp89;
  const Scalar _tmp131 = _tmp49 + _tmp83;
  const Scalar _tmp132 = _tmp131 * _tmp86;
  const Scalar _tmp133 = Scalar(1.0) / (-_tmp132 - _tmp45 + _tmp85);
  const Scalar _tmp134 = _tmp131 * _tmp133;
  const Scalar _tmp135 = _tmp101 * _tmp134 - _tmp102 * _tmp130 - _tmp95;
  const Scalar _tmp136 = Scalar(1.0) * _tmp84;
  const Scalar _tmp137 = Scalar(1.0) * _tmp133;
  const Scalar _tmp138 = Scalar(1.0) * _tmp109 * (_tmp101 * _tmp137 - _tmp135 * _tmp136);
  const Scalar _tmp139 = _tmp54 + _tmp78;
  const Scalar _tmp140 = _tmp139 * fh1;
  const Scalar _tmp141 = -_tmp109 * _tmp140 - Scalar(40.024799999999999) * _tmp11 - _tmp17 * fv1;
  const Scalar _tmp142 = _tmp84 * (_tmp132 * _tmp137 + Scalar(1.0));
  const Scalar _tmp143 = _tmp137 * _tmp86;
  const Scalar _tmp144 = -Scalar(1.0) * _tmp142 + Scalar(1.0) * _tmp143;
  const Scalar _tmp145 = _tmp116 * _tmp140 + Scalar(40.024799999999999) * _tmp23 + _tmp27 * fv1;
  const Scalar _tmp146 = _tmp131 * _tmp137;
  const Scalar _tmp147 = _tmp146 * _tmp84;
  const Scalar _tmp148 = -Scalar(1.0) * _tmp137 + Scalar(1.0) * _tmp147;
  const Scalar _tmp149 = _tmp121 + _tmp123 * _tmp134 - _tmp124 * _tmp130;
  const Scalar _tmp150 = Scalar(1.0) * _tmp116 * (_tmp123 * _tmp137 - _tmp136 * _tmp149);
  const Scalar _tmp151 = _tmp100 * _tmp146 - _tmp118 * _tmp81;
  const Scalar _tmp152 = _tmp100 * _tmp137;
  const Scalar _tmp153 = Scalar(1.0) * _tmp117;
  const Scalar _tmp154 = _tmp153 * (-_tmp136 * _tmp151 + _tmp152);
  const Scalar _tmp155 = _tmp112 * _tmp133;
  const Scalar _tmp156 = _tmp84 * (-_tmp113 * _tmp130 - _tmp131 * _tmp155 + _tmp83);
  const Scalar _tmp157 =
      Scalar(1.0) * _tmp111 * (-_tmp112 * _tmp137 - Scalar(1.0) * _tmp156 + Scalar(1.0)) +
      _tmp138 * fh1 + _tmp141 * _tmp144 + _tmp145 * _tmp148 + _tmp150 * fh1 + _tmp154 * fh1;
  const Scalar _tmp158 = std::asinh(_tmp129 * _tmp157);
  const Scalar _tmp159 = Scalar(1.4083112389913199) * _tmp128;
  const Scalar _tmp160 =
      -_tmp158 * _tmp159 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp41), Scalar(2)) +
                     Scalar(0.77439999999999998) *
                         std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp43), Scalar(2))));
  const Scalar _tmp161 = Scalar(0.71007031138673404) * _tmp129;
  const Scalar _tmp162 = _tmp160 * _tmp161;
  const Scalar _tmp163 = Scalar(1.0) * _tmp158;
  const Scalar _tmp164 = -_tmp110 - _tmp120 - _tmp127;
  const Scalar _tmp165 = Scalar(1.4083112389913199) * _tmp164;
  const Scalar _tmp166 = std::pow(_tmp128, Scalar(-2));
  const Scalar _tmp167 = _tmp164 * _tmp166;
  const Scalar _tmp168 = _tmp109 * _tmp139;
  const Scalar _tmp169 = _tmp116 * _tmp139;
  const Scalar _tmp170 =
      (_tmp129 * (_tmp138 - _tmp144 * _tmp168 + _tmp148 * _tmp169 + _tmp150 + _tmp154) -
       _tmp157 * _tmp167) /
      std::sqrt(Scalar(std::pow(_tmp157, Scalar(2)) * _tmp166 + 1));
  const Scalar _tmp171 = _tmp109 * _tmp135 * _tmp84;
  const Scalar _tmp172 = _tmp137 * _tmp145;
  const Scalar _tmp173 = _tmp116 * _tmp149 * _tmp84;
  const Scalar _tmp174 = _tmp117 * _tmp151 * _tmp84;
  const Scalar _tmp175 = _tmp111 * _tmp156 - _tmp131 * _tmp172 * _tmp84 + _tmp141 * _tmp142 +
                         _tmp171 * fh1 + _tmp173 * fh1 + _tmp174 * fh1;
  const Scalar _tmp176 = _tmp153 * _tmp99;
  const Scalar _tmp177 = _tmp176 * fh1;
  const Scalar _tmp178 = _tmp105 * _tmp109;
  const Scalar _tmp179 = _tmp116 * _tmp126;
  const Scalar _tmp180 = -_tmp111 * _tmp115 - _tmp119 * _tmp177 + _tmp178 * fh1 + _tmp179 * fh1;
  const Scalar _tmp181 = Scalar(1.0) / (_tmp180);
  const Scalar _tmp182 = std::asinh(_tmp175 * _tmp181);
  const Scalar _tmp183 = Scalar(1.0) * _tmp182;
  const Scalar _tmp184 = Scalar(1.4083112389913199) * _tmp180;
  const Scalar _tmp185 =
      -_tmp182 * _tmp184 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(Scalar(0.77439999999999998) *
                         std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp69 - 1), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp72 - 1), Scalar(2))));
  const Scalar _tmp186 = Scalar(0.71007031138673404) * _tmp181;
  const Scalar _tmp187 = _tmp185 * _tmp186;
  const Scalar _tmp188 = -_tmp119 * _tmp176 + _tmp178 + _tmp179;
  const Scalar _tmp189 = Scalar(1.4083112389913199) * _tmp188;
  const Scalar _tmp190 = std::pow(_tmp180, Scalar(-2));
  const Scalar _tmp191 = _tmp188 * _tmp190;
  const Scalar _tmp192 = (-_tmp175 * _tmp191 + _tmp181 * (-_tmp142 * _tmp168 - _tmp147 * _tmp169 +
                                                          _tmp171 + _tmp173 + _tmp174)) /
                         std::sqrt(Scalar(std::pow(_tmp175, Scalar(2)) * _tmp190 + 1));
  const Scalar _tmp193 = _tmp116 * _tmp123 * _tmp133;
  const Scalar _tmp194 = _tmp117 * _tmp152;
  const Scalar _tmp195 = _tmp101 * _tmp109 * _tmp133;
  const Scalar _tmp196 = _tmp111 * _tmp155 - _tmp141 * _tmp143 + _tmp172 - _tmp193 * fh1 -
                         _tmp194 * fh1 - _tmp195 * fh1;
  const Scalar _tmp197 = _tmp116 * _tmp125;
  const Scalar _tmp198 = _tmp103 * _tmp109;
  const Scalar _tmp199 = _tmp111 * _tmp114 + _tmp177 + _tmp197 * fh1 + _tmp198 * fh1;
  const Scalar _tmp200 = Scalar(1.0) / (_tmp199);
  const Scalar _tmp201 = std::asinh(_tmp196 * _tmp200);
  const Scalar _tmp202 = Scalar(1.0) * _tmp201;
  const Scalar _tmp203 = Scalar(1.4083112389913199) * _tmp199;
  const Scalar _tmp204 =
      -_tmp201 * _tmp203 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(Scalar(0.77439999999999998) *
                         std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp50), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp46 - 1), Scalar(2))));
  const Scalar _tmp205 = Scalar(0.71007031138673404) * _tmp200;
  const Scalar _tmp206 = _tmp204 * _tmp205;
  const Scalar _tmp207 = _tmp176 + _tmp197 + _tmp198;
  const Scalar _tmp208 = Scalar(1.4083112389913199) * _tmp207;
  const Scalar _tmp209 = std::pow(_tmp199, Scalar(-2));
  const Scalar _tmp210 = _tmp207 * _tmp209;
  const Scalar _tmp211 = (-_tmp196 * _tmp210 + _tmp200 * (_tmp137 * _tmp169 + _tmp143 * _tmp168 -
                                                          _tmp193 - _tmp194 - _tmp195)) /
                         std::sqrt(Scalar(std::pow(_tmp196, Scalar(2)) * _tmp209 + 1));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = Scalar(1.4083112389913199) * fh1 *
                   (Scalar(1.0) * _tmp37 * _tmp38 * fv1 * std::cosh(_tmp36) -
                    (-Scalar(0.71007031138673404) * _tmp33 * _tmp37 +
                     _tmp34 * (Scalar(1.4083112389913199) * _tmp30 * _tmp38 - _tmp32)) *
                        std::cosh(_tmp35)) -
               Scalar(1.4083112389913199) * std::sinh(_tmp35) -
               Scalar(1.4083112389913199) * std::sinh(_tmp36);
  _res(1, 0) = _tmp159 * (-Scalar(1.0) * _tmp170 * std::cosh(_tmp163) -
                          (-Scalar(0.71007031138673404) * _tmp160 * _tmp167 +
                           _tmp161 * (-_tmp158 * _tmp165 - _tmp159 * _tmp170)) *
                              std::cosh(_tmp162)) +
               _tmp165 * (-std::sinh(_tmp162) - std::sinh(_tmp163));
  _res(2, 0) = _tmp184 * (-Scalar(1.0) * _tmp192 * std::cosh(_tmp183) -
                          (-Scalar(0.71007031138673404) * _tmp185 * _tmp191 +
                           _tmp186 * (-_tmp182 * _tmp189 - _tmp184 * _tmp192)) *
                              std::cosh(_tmp187)) +
               _tmp189 * (-std::sinh(_tmp183) - std::sinh(_tmp187));
  _res(3, 0) = _tmp203 * (-Scalar(1.0) * _tmp211 * std::cosh(_tmp202) -
                          (-Scalar(0.71007031138673404) * _tmp204 * _tmp210 +
                           _tmp205 * (-_tmp201 * _tmp208 - _tmp203 * _tmp211)) *
                              std::cosh(_tmp206)) +
               _tmp208 * (-std::sinh(_tmp202) - std::sinh(_tmp206));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
