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
 * Symbolic function: IK_residual_func_cost2_wrt_fv1_Nl12
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2WrtFv1Nl12(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 597

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (191)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp2 = Scalar(1.0) * _tmp0 /
                       std::sqrt(Scalar(1 + std::pow(fv1, Scalar(2)) / std::pow(fh1, Scalar(2))));
  const Scalar _tmp3 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp4 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp6 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp7 = Scalar(0.20999999999999999) * _tmp4 + Scalar(0.20999999999999999) * _tmp6 +
                       Scalar(0.20999999999999999);
  const Scalar _tmp8 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp9 = 2 * _tmp5 * _tmp8;
  const Scalar _tmp10 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                        2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp11 = _tmp10 * _tmp3;
  const Scalar _tmp12 = _tmp11 + _tmp9;
  const Scalar _tmp13 = -Scalar(0.010999999999999999) * _tmp12;
  const Scalar _tmp14 = 2 * _tmp3;
  const Scalar _tmp15 = _tmp14 * _tmp8;
  const Scalar _tmp16 = _tmp10 * _tmp5;
  const Scalar _tmp17 = Scalar(0.20999999999999999) * _tmp15 - Scalar(0.20999999999999999) * _tmp16;
  const Scalar _tmp18 = _tmp13 + _tmp17;
  const Scalar _tmp19 = _tmp18 + _tmp7;
  const Scalar _tmp20 = _tmp19 + position_vector(0, 0);
  const Scalar _tmp21 = 1 - 2 * std::pow(_tmp8, Scalar(2));
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp21 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp23 = _tmp14 * _tmp5;
  const Scalar _tmp24 = _tmp10 * _tmp8;
  const Scalar _tmp25 = _tmp23 - _tmp24;
  const Scalar _tmp26 = Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = -_tmp26;
  const Scalar _tmp28 = Scalar(0.20999999999999999) * _tmp15 + Scalar(0.20999999999999999) * _tmp16;
  const Scalar _tmp29 = _tmp27 + _tmp28;
  const Scalar _tmp30 = _tmp22 + _tmp29;
  const Scalar _tmp31 = _tmp30 + position_vector(1, 0);
  const Scalar _tmp32 = Scalar(1.4083112389913199) * fh1;
  const Scalar _tmp33 = -_tmp22;
  const Scalar _tmp34 = -_tmp28;
  const Scalar _tmp35 = _tmp27 + _tmp34;
  const Scalar _tmp36 = _tmp33 + _tmp35;
  const Scalar _tmp37 = _tmp36 + position_vector(1, 0);
  const Scalar _tmp38 = -_tmp7;
  const Scalar _tmp39 = _tmp13 - _tmp17;
  const Scalar _tmp40 = _tmp38 + _tmp39;
  const Scalar _tmp41 = _tmp40 + position_vector(0, 0);
  const Scalar _tmp42 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp43 = -Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp44 = -_tmp43;
  const Scalar _tmp45 =
      -Scalar(0.010999999999999999) * _tmp21 - Scalar(0.010999999999999999) * _tmp4;
  const Scalar _tmp46 = Scalar(0.20999999999999999) * _tmp23 + Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp47 = _tmp45 + _tmp46;
  const Scalar _tmp48 = _tmp44 + _tmp47;
  const Scalar _tmp49 = _tmp18 + _tmp38;
  const Scalar _tmp50 = _tmp49 + position_vector(0, 0);
  const Scalar _tmp51 = _tmp50 + Scalar(125.0);
  const Scalar _tmp52 = _tmp22 + _tmp35;
  const Scalar _tmp53 = _tmp52 + position_vector(1, 0);
  const Scalar _tmp54 = _tmp53 + Scalar(-110.0);
  const Scalar _tmp55 = std::pow(Scalar(std::pow(_tmp51, Scalar(2)) + std::pow(_tmp54, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp56 = _tmp51 * _tmp55;
  const Scalar _tmp57 = _tmp37 + Scalar(110.0);
  const Scalar _tmp58 = _tmp41 + Scalar(125.0);
  const Scalar _tmp59 = Scalar(1.0) / (_tmp58);
  const Scalar _tmp60 = _tmp57 * _tmp59;
  const Scalar _tmp61 = _tmp45 - _tmp46;
  const Scalar _tmp62 = _tmp44 + _tmp61;
  const Scalar _tmp63 = _tmp56 * _tmp62;
  const Scalar _tmp64 = _tmp54 * _tmp55;
  const Scalar _tmp65 = _tmp43 + _tmp61;
  const Scalar _tmp66 = _tmp29 + _tmp33;
  const Scalar _tmp67 = _tmp66 + position_vector(1, 0);
  const Scalar _tmp68 = _tmp67 + Scalar(110.0);
  const Scalar _tmp69 = _tmp39 + _tmp7;
  const Scalar _tmp70 = _tmp69 + position_vector(0, 0);
  const Scalar _tmp71 = _tmp70 + Scalar(-125.0);
  const Scalar _tmp72 = std::pow(Scalar(std::pow(_tmp68, Scalar(2)) + std::pow(_tmp71, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp73 = _tmp68 * _tmp72;
  const Scalar _tmp74 = _tmp71 * _tmp72;
  const Scalar _tmp75 = _tmp60 * _tmp62;
  const Scalar _tmp76 = _tmp65 * _tmp73 - _tmp74 * _tmp75;
  const Scalar _tmp77 = Scalar(1.0) / (_tmp60 * _tmp74 - _tmp73);
  const Scalar _tmp78 = _tmp56 * _tmp60 - _tmp64;
  const Scalar _tmp79 = _tmp77 * _tmp78;
  const Scalar _tmp80 = _tmp48 * _tmp64 - _tmp60 * _tmp63 - _tmp76 * _tmp79;
  const Scalar _tmp81 = Scalar(1.0) * _tmp36;
  const Scalar _tmp82 = -_tmp81;
  const Scalar _tmp83 = Scalar(1.0) / (_tmp66 + _tmp82);
  const Scalar _tmp84 = Scalar(1.0) * _tmp40;
  const Scalar _tmp85 = _tmp83 * (-_tmp69 + _tmp84);
  const Scalar _tmp86 = _tmp62 * _tmp74 - _tmp65 * _tmp74;
  const Scalar _tmp87 = -_tmp48 * _tmp56 + _tmp63 - _tmp79 * _tmp86 - _tmp80 * _tmp85;
  const Scalar _tmp88 = Scalar(1.0) / (_tmp87);
  const Scalar _tmp89 = _tmp81 * _tmp85 + _tmp84;
  const Scalar _tmp90 = 0;
  const Scalar _tmp91 = _tmp79 * _tmp90;
  const Scalar _tmp92 =
      std::sqrt(Scalar(std::pow(_tmp57, Scalar(2)) + std::pow(_tmp58, Scalar(2))));
  const Scalar _tmp93 = _tmp59 * _tmp92;
  const Scalar _tmp94 = _tmp93 * (_tmp56 * _tmp90 - _tmp74 * _tmp91);
  const Scalar _tmp95 = Scalar(1.0) / (_tmp92);
  const Scalar _tmp96 = _tmp93 * (-_tmp36 * _tmp58 * _tmp95 + _tmp40 * _tmp57 * _tmp95);
  const Scalar _tmp97 = _tmp66 * _tmp74 - _tmp69 * _tmp73 + _tmp74 * _tmp96;
  const Scalar _tmp98 = -_tmp49 * _tmp64 + _tmp52 * _tmp56 + _tmp56 * _tmp96 - _tmp79 * _tmp97;
  const Scalar _tmp99 = Scalar(1.0) / (_tmp98);
  const Scalar _tmp100 = Scalar(1.0) * _tmp77;
  const Scalar _tmp101 = _tmp100 * _tmp78 * _tmp99;
  const Scalar _tmp102 = Scalar(1.0) * _tmp99;
  const Scalar _tmp103 = _tmp31 + Scalar(-110.0);
  const Scalar _tmp104 = _tmp20 + Scalar(-125.0);
  const Scalar _tmp105 =
      std::pow(Scalar(std::pow(_tmp103, Scalar(2)) + std::pow(_tmp104, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp106 = _tmp104 * _tmp105;
  const Scalar _tmp107 = _tmp103 * _tmp105;
  const Scalar _tmp108 = fh1 * (-_tmp106 * _tmp30 + _tmp107 * _tmp19);
  const Scalar _tmp109 = _tmp100 * _tmp76;
  const Scalar _tmp110 = -_tmp100 * _tmp86 + _tmp109 * _tmp85;
  const Scalar _tmp111 = _tmp88 * _tmp98;
  const Scalar _tmp112 = _tmp87 * _tmp99;
  const Scalar _tmp113 = _tmp112 * (-_tmp100 * _tmp97 - _tmp110 * _tmp111);
  const Scalar _tmp114 = _tmp88 * (_tmp110 + _tmp113);
  const Scalar _tmp115 = -_tmp114 * _tmp78 + Scalar(1.0);
  const Scalar _tmp116 = _tmp74 * _tmp77;
  const Scalar _tmp117 = _tmp107 * fh1;
  const Scalar _tmp118 = _tmp60 * _tmp77;
  const Scalar _tmp119 = _tmp118 * _tmp76 + _tmp75;
  const Scalar _tmp120 = _tmp118 * _tmp86 - _tmp119 * _tmp85 - _tmp62;
  const Scalar _tmp121 = _tmp112 * (-_tmp111 * _tmp120 + _tmp118 * _tmp97 - _tmp96);
  const Scalar _tmp122 = _tmp88 * (_tmp120 + _tmp121);
  const Scalar _tmp123 = -_tmp122 * _tmp78 - _tmp60;
  const Scalar _tmp124 = _tmp106 * fh1;
  const Scalar _tmp125 = -_tmp108 * _tmp93 * (-_tmp101 * _tmp74 + _tmp102 * _tmp56) -
                         _tmp117 * _tmp93 * (_tmp114 * _tmp56 + _tmp115 * _tmp116) -
                         _tmp124 * _tmp93 * (_tmp116 * _tmp123 + _tmp122 * _tmp56 + Scalar(1.0)) -
                         _tmp42 * _tmp94;
  const Scalar _tmp126 = Scalar(1.0) / (_tmp125);
  const Scalar _tmp127 = _tmp52 + _tmp82;
  const Scalar _tmp128 = _tmp127 * _tmp85;
  const Scalar _tmp129 = Scalar(1.0) / (-_tmp128 - _tmp49 + _tmp84);
  const Scalar _tmp130 = _tmp127 * _tmp129;
  const Scalar _tmp131 = -_tmp109 + _tmp113 * _tmp130 - _tmp114 * _tmp80;
  const Scalar _tmp132 = Scalar(1.0) * _tmp83;
  const Scalar _tmp133 = Scalar(1.0) * _tmp129;
  const Scalar _tmp134 = _tmp112 * _tmp133;
  const Scalar _tmp135 = -_tmp102 * _tmp80 + _tmp127 * _tmp134;
  const Scalar _tmp136 = _tmp129 * _tmp89;
  const Scalar _tmp137 = _tmp83 * (-_tmp127 * _tmp136 - _tmp80 * _tmp90 + _tmp82);
  const Scalar _tmp138 = -Scalar(1.0) * _tmp133 * _tmp89 - Scalar(1.0) * _tmp137 + Scalar(1.0);
  const Scalar _tmp139 = fh1 * (_tmp43 + _tmp47);
  const Scalar _tmp140 = -_tmp107 * _tmp139 - Scalar(40.024799999999999) * _tmp25 - _tmp30 * fv1;
  const Scalar _tmp141 = _tmp133 * _tmp85;
  const Scalar _tmp142 = _tmp128 * _tmp133 + Scalar(1.0);
  const Scalar _tmp143 = -Scalar(1.0) * _tmp132 * _tmp142 + Scalar(1.0) * _tmp141;
  const Scalar _tmp144 = _tmp119 + _tmp121 * _tmp130 - _tmp122 * _tmp80;
  const Scalar _tmp145 = _tmp106 * _tmp139 + Scalar(40.024799999999999) * _tmp12 + _tmp19 * fv1;
  const Scalar _tmp146 = _tmp127 * _tmp83;
  const Scalar _tmp147 = Scalar(1.0) * _tmp133 * _tmp146 - Scalar(1.0) * _tmp133;
  const Scalar _tmp148 = Scalar(1.0) * _tmp108 * (-_tmp132 * _tmp135 + _tmp134) +
                         Scalar(1.0) * _tmp117 * (_tmp113 * _tmp133 - _tmp131 * _tmp132) +
                         Scalar(1.0) * _tmp124 * (_tmp121 * _tmp133 - _tmp132 * _tmp144) +
                         _tmp138 * _tmp42 + _tmp140 * _tmp143 + _tmp145 * _tmp147;
  const Scalar _tmp149 = std::asinh(_tmp126 * _tmp148);
  const Scalar _tmp150 = Scalar(1.4083112389913199) * _tmp125;
  const Scalar _tmp151 =
      -_tmp149 * _tmp150 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(Scalar(0.77439999999999998) *
                         std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp37 - 1), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp41 - 1), Scalar(2))));
  const Scalar _tmp152 = std::pow(_tmp125, Scalar(-2));
  const Scalar _tmp153 = _tmp152 * _tmp94;
  const Scalar _tmp154 = Scalar(1.4083112389913199) * _tmp94;
  const Scalar _tmp155 = _tmp26 + _tmp33 + _tmp34;
  const Scalar _tmp156 =
      (_tmp126 * (-_tmp138 + _tmp143 * _tmp155 + _tmp147 * _tmp19) - _tmp148 * _tmp153) /
      std::sqrt(Scalar(std::pow(_tmp148, Scalar(2)) * _tmp152 + 1));
  const Scalar _tmp157 = Scalar(0.71007031138673404) * _tmp126;
  const Scalar _tmp158 = _tmp151 * _tmp157;
  const Scalar _tmp159 = Scalar(1.0) * _tmp149;
  const Scalar _tmp160 = _tmp42 * _tmp90;
  const Scalar _tmp161 = -_tmp101 * _tmp108 + _tmp115 * _tmp117 * _tmp77 +
                         _tmp123 * _tmp124 * _tmp77 - _tmp160 * _tmp79;
  const Scalar _tmp162 = Scalar(1.0) / (_tmp161);
  const Scalar _tmp163 = _tmp83 * fh1;
  const Scalar _tmp164 = _tmp133 * _tmp145;
  const Scalar _tmp165 = _tmp142 * _tmp83;
  const Scalar _tmp166 = _tmp106 * _tmp144 * _tmp163 + _tmp107 * _tmp131 * _tmp163 +
                         _tmp108 * _tmp135 * _tmp83 + _tmp137 * _tmp42 + _tmp140 * _tmp165 -
                         _tmp146 * _tmp164;
  const Scalar _tmp167 = std::asinh(_tmp162 * _tmp166);
  const Scalar _tmp168 = Scalar(1.4083112389913199) * _tmp161;
  const Scalar _tmp169 =
      -_tmp167 * _tmp168 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp70), Scalar(2)) +
                     Scalar(0.77439999999999998) *
                         std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp67 - 1), Scalar(2))));
  const Scalar _tmp170 = Scalar(0.71007031138673404) * _tmp162;
  const Scalar _tmp171 = _tmp169 * _tmp170;
  const Scalar _tmp172 = Scalar(1.0) * _tmp167;
  const Scalar _tmp173 = Scalar(1.4083112389913199) * _tmp90;
  const Scalar _tmp174 = _tmp173 * _tmp79;
  const Scalar _tmp175 = std::pow(_tmp161, Scalar(-2));
  const Scalar _tmp176 = _tmp175 * _tmp91;
  const Scalar _tmp177 = _tmp133 * _tmp19;
  const Scalar _tmp178 =
      (_tmp162 * (-_tmp137 - _tmp146 * _tmp177 + _tmp155 * _tmp165) - _tmp166 * _tmp176) /
      std::sqrt(Scalar(std::pow(_tmp166, Scalar(2)) * _tmp175 + 1));
  const Scalar _tmp179 = -_tmp108 * _tmp134 - _tmp113 * _tmp117 * _tmp129 -
                         _tmp121 * _tmp124 * _tmp129 + _tmp136 * _tmp42 - _tmp140 * _tmp141 +
                         _tmp164;
  const Scalar _tmp180 = _tmp102 * _tmp108 + _tmp114 * _tmp117 + _tmp122 * _tmp124 + _tmp160;
  const Scalar _tmp181 = Scalar(1.0) / (_tmp180);
  const Scalar _tmp182 = std::asinh(_tmp179 * _tmp181);
  const Scalar _tmp183 = Scalar(1.4083112389913199) * _tmp180;
  const Scalar _tmp184 =
      -_tmp182 * _tmp183 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(Scalar(0.77439999999999998) *
                         std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp53), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp50 - 1), Scalar(2))));
  const Scalar _tmp185 = Scalar(0.71007031138673404) * _tmp181;
  const Scalar _tmp186 = _tmp184 * _tmp185;
  const Scalar _tmp187 = Scalar(1.0) * _tmp182;
  const Scalar _tmp188 = std::pow(_tmp180, Scalar(-2));
  const Scalar _tmp189 = _tmp188 * _tmp90;
  const Scalar _tmp190 = (_tmp179 * _tmp189 + _tmp181 * (-_tmp136 - _tmp141 * _tmp155 + _tmp177)) /
                         std::sqrt(Scalar(std::pow(_tmp179, Scalar(2)) * _tmp188 + 1));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      _tmp32 *
      (-_tmp2 * std::cosh(Scalar(1.0) * _tmp1) +
       _tmp2 * std::cosh(Scalar(0.71007031138673404) * _tmp0 *
                         (-_tmp1 * _tmp32 -
                          Scalar(125.0) *
                              std::sqrt(Scalar(
                                  std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp20),
                                           Scalar(2)) +
                                  Scalar(0.77439999999999998) *
                                      std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp31),
                                               Scalar(2)))))));
  _res(1, 0) = _tmp150 * (-Scalar(1.0) * _tmp156 * std::cosh(_tmp159) -
                          (-Scalar(0.71007031138673404) * _tmp151 * _tmp153 +
                           _tmp157 * (-_tmp149 * _tmp154 - _tmp150 * _tmp156)) *
                              std::cosh(_tmp158)) +
               _tmp154 * (-std::sinh(_tmp158) - std::sinh(_tmp159));
  _res(2, 0) = _tmp168 * (-Scalar(1.0) * _tmp178 * std::cosh(_tmp172) -
                          (-Scalar(0.71007031138673404) * _tmp169 * _tmp176 +
                           _tmp170 * (-_tmp167 * _tmp174 - _tmp168 * _tmp178)) *
                              std::cosh(_tmp171)) +
               _tmp174 * (-std::sinh(_tmp171) - std::sinh(_tmp172));
  _res(3, 0) = -_tmp173 * (-std::sinh(_tmp186) - std::sinh(_tmp187)) +
               _tmp183 * (-Scalar(1.0) * _tmp190 * std::cosh(_tmp187) -
                          (Scalar(0.71007031138673404) * _tmp184 * _tmp189 +
                           _tmp185 * (_tmp173 * _tmp182 - _tmp183 * _tmp190)) *
                              std::cosh(_tmp186));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
