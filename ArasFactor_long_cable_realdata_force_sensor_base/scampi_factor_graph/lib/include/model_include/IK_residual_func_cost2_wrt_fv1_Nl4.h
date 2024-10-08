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
 * Symbolic function: IK_residual_func_cost2_wrt_fv1_Nl4
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2WrtFv1Nl4(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 590

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (189)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp2 = Scalar(1.0) * _tmp0 /
                       std::sqrt(Scalar(1 + std::pow(fv1, Scalar(2)) / std::pow(fh1, Scalar(2))));
  const Scalar _tmp3 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp4 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp6 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp7 = Scalar(0.20999999999999999) * _tmp4 + Scalar(0.20999999999999999) * _tmp6 +
                       Scalar(0.20999999999999999);
  const Scalar _tmp8 = -_tmp7;
  const Scalar _tmp9 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp10 = 2 * _tmp3;
  const Scalar _tmp11 = _tmp10 * _tmp9;
  const Scalar _tmp12 = -_DeltaRot[0] * _Rot_init[0] - _DeltaRot[1] * _Rot_init[1] -
                        _DeltaRot[2] * _Rot_init[2] + _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp13 = 2 * _tmp12 * _tmp5;
  const Scalar _tmp14 = Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp15 = 2 * _tmp9;
  const Scalar _tmp16 = _tmp15 * _tmp5;
  const Scalar _tmp17 = _tmp10 * _tmp12;
  const Scalar _tmp18 = _tmp16 - _tmp17;
  const Scalar _tmp19 = Scalar(0.010999999999999999) * _tmp18;
  const Scalar _tmp20 = -_tmp19;
  const Scalar _tmp21 = -_tmp14 + _tmp20;
  const Scalar _tmp22 = _tmp21 + _tmp8;
  const Scalar _tmp23 = _tmp22 + position_vector(1, 0);
  const Scalar _tmp24 = 1 - 2 * std::pow(_tmp9, Scalar(2));
  const Scalar _tmp25 = Scalar(0.20999999999999999) * _tmp24 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp26 = -_tmp25;
  const Scalar _tmp27 = _tmp10 * _tmp5;
  const Scalar _tmp28 = _tmp12 * _tmp15;
  const Scalar _tmp29 = _tmp27 + _tmp28;
  const Scalar _tmp30 = -Scalar(0.010999999999999999) * _tmp29;
  const Scalar _tmp31 = Scalar(0.20999999999999999) * _tmp11 - Scalar(0.20999999999999999) * _tmp13;
  const Scalar _tmp32 = _tmp30 - _tmp31;
  const Scalar _tmp33 = _tmp26 + _tmp32;
  const Scalar _tmp34 = _tmp33 + position_vector(0, 0);
  const Scalar _tmp35 = Scalar(1.4083112389913199) * fh1;
  const Scalar _tmp36 = _tmp14 + _tmp7;
  const Scalar _tmp37 = _tmp20 + _tmp36;
  const Scalar _tmp38 = _tmp37 + position_vector(1, 0);
  const Scalar _tmp39 = _tmp38 + Scalar(-110.0);
  const Scalar _tmp40 = _tmp30 + _tmp31;
  const Scalar _tmp41 = _tmp25 + _tmp40;
  const Scalar _tmp42 = _tmp41 + position_vector(0, 0);
  const Scalar _tmp43 = _tmp42 + Scalar(-125.0);
  const Scalar _tmp44 = std::pow(Scalar(std::pow(_tmp39, Scalar(2)) + std::pow(_tmp43, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp45 = _tmp43 * _tmp44;
  const Scalar _tmp46 = _tmp39 * _tmp44;
  const Scalar _tmp47 = _tmp14 + _tmp20 + _tmp8;
  const Scalar _tmp48 = _tmp47 + position_vector(1, 0);
  const Scalar _tmp49 = _tmp48 + Scalar(110.0);
  const Scalar _tmp50 = _tmp25 + _tmp32;
  const Scalar _tmp51 = _tmp50 + position_vector(0, 0);
  const Scalar _tmp52 = _tmp51 + Scalar(-125.0);
  const Scalar _tmp53 = std::pow(Scalar(std::pow(_tmp49, Scalar(2)) + std::pow(_tmp52, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp54 = _tmp52 * _tmp53;
  const Scalar _tmp55 = _tmp21 + _tmp7;
  const Scalar _tmp56 = _tmp55 + position_vector(1, 0);
  const Scalar _tmp57 = _tmp56 + Scalar(-110.0);
  const Scalar _tmp58 = _tmp26 + _tmp40;
  const Scalar _tmp59 = _tmp58 + position_vector(0, 0);
  const Scalar _tmp60 = _tmp59 + Scalar(125.0);
  const Scalar _tmp61 =
      std::sqrt(Scalar(std::pow(_tmp57, Scalar(2)) + std::pow(_tmp60, Scalar(2))));
  const Scalar _tmp62 = Scalar(1.0) / (_tmp61);
  const Scalar _tmp63 = Scalar(1.0) / (_tmp60);
  const Scalar _tmp64 = _tmp61 * _tmp63;
  const Scalar _tmp65 = _tmp64 * (-_tmp55 * _tmp60 * _tmp62 + _tmp57 * _tmp58 * _tmp62);
  const Scalar _tmp66 = _tmp49 * _tmp53;
  const Scalar _tmp67 = _tmp47 * _tmp54 - _tmp50 * _tmp66 + _tmp54 * _tmp65;
  const Scalar _tmp68 = _tmp57 * _tmp63;
  const Scalar _tmp69 = Scalar(1.0) / (_tmp54 * _tmp68 - _tmp66);
  const Scalar _tmp70 = _tmp45 * _tmp68 - _tmp46;
  const Scalar _tmp71 = _tmp69 * _tmp70;
  const Scalar _tmp72 = _tmp37 * _tmp45 - _tmp41 * _tmp46 + _tmp45 * _tmp65 - _tmp67 * _tmp71;
  const Scalar _tmp73 = Scalar(1.0) / (_tmp72);
  const Scalar _tmp74 = Scalar(1.0) * _tmp73;
  const Scalar _tmp75 = Scalar(1.0) * _tmp69;
  const Scalar _tmp76 = _tmp70 * _tmp73 * _tmp75;
  const Scalar _tmp77 = _tmp23 + Scalar(110.0);
  const Scalar _tmp78 = _tmp34 + Scalar(125.0);
  const Scalar _tmp79 = std::pow(Scalar(std::pow(_tmp77, Scalar(2)) + std::pow(_tmp78, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp80 = _tmp77 * _tmp79;
  const Scalar _tmp81 = _tmp78 * _tmp79;
  const Scalar _tmp82 = fh1 * (-_tmp22 * _tmp81 + _tmp33 * _tmp80);
  const Scalar _tmp83 = Scalar(0.20999999999999999) * _tmp16 + Scalar(0.20999999999999999) * _tmp17;
  const Scalar _tmp84 =
      -Scalar(0.010999999999999999) * _tmp24 - Scalar(0.010999999999999999) * _tmp4;
  const Scalar _tmp85 = Scalar(0.20999999999999999) * _tmp27 - Scalar(0.20999999999999999) * _tmp28;
  const Scalar _tmp86 = _tmp84 + _tmp85;
  const Scalar _tmp87 = _tmp83 + _tmp86;
  const Scalar _tmp88 = _tmp84 - _tmp85;
  const Scalar _tmp89 = _tmp83 + _tmp88;
  const Scalar _tmp90 = _tmp45 * _tmp89;
  const Scalar _tmp91 = -_tmp83;
  const Scalar _tmp92 = _tmp86 + _tmp91;
  const Scalar _tmp93 = _tmp54 * _tmp89;
  const Scalar _tmp94 = _tmp66 * _tmp92 - _tmp68 * _tmp93;
  const Scalar _tmp95 = _tmp46 * _tmp87 - _tmp68 * _tmp90 - _tmp71 * _tmp94;
  const Scalar _tmp96 = Scalar(1.0) * _tmp55;
  const Scalar _tmp97 = -_tmp96;
  const Scalar _tmp98 = Scalar(1.0) / (_tmp47 + _tmp97);
  const Scalar _tmp99 = Scalar(1.0) * _tmp58;
  const Scalar _tmp100 = _tmp98 * (-_tmp50 + _tmp99);
  const Scalar _tmp101 = -_tmp54 * _tmp92 + _tmp93;
  const Scalar _tmp102 = -_tmp100 * _tmp95 - _tmp101 * _tmp71 - _tmp45 * _tmp87 + _tmp90;
  const Scalar _tmp103 = Scalar(1.0) / (_tmp102);
  const Scalar _tmp104 = _tmp75 * _tmp94;
  const Scalar _tmp105 = _tmp100 * _tmp104 - _tmp101 * _tmp75;
  const Scalar _tmp106 = _tmp103 * _tmp72;
  const Scalar _tmp107 = _tmp102 * _tmp73;
  const Scalar _tmp108 = _tmp107 * (-_tmp105 * _tmp106 - _tmp67 * _tmp75);
  const Scalar _tmp109 = _tmp103 * (_tmp105 + _tmp108);
  const Scalar _tmp110 = -_tmp109 * _tmp70 + Scalar(1.0);
  const Scalar _tmp111 = _tmp54 * _tmp69;
  const Scalar _tmp112 = _tmp80 * fh1;
  const Scalar _tmp113 = _tmp68 * _tmp69;
  const Scalar _tmp114 = _tmp113 * _tmp94 + _tmp68 * _tmp89;
  const Scalar _tmp115 = -_tmp100 * _tmp114 + _tmp101 * _tmp113 - _tmp89;
  const Scalar _tmp116 = _tmp107 * (-_tmp106 * _tmp115 + _tmp113 * _tmp67 - _tmp65);
  const Scalar _tmp117 = _tmp103 * (_tmp115 + _tmp116);
  const Scalar _tmp118 = -_tmp117 * _tmp70 - _tmp68;
  const Scalar _tmp119 = _tmp81 * fh1;
  const Scalar _tmp120 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp121 = _tmp100 * _tmp96 + _tmp99;
  const Scalar _tmp122 = 0;
  const Scalar _tmp123 = _tmp122 * _tmp71;
  const Scalar _tmp124 = _tmp64 * (_tmp122 * _tmp45 - _tmp123 * _tmp54);
  const Scalar _tmp125 = -_tmp112 * _tmp64 * (_tmp109 * _tmp45 + _tmp110 * _tmp111) -
                         _tmp119 * _tmp64 * (_tmp111 * _tmp118 + _tmp117 * _tmp45 + Scalar(1.0)) -
                         _tmp120 * _tmp124 - _tmp64 * _tmp82 * (_tmp45 * _tmp74 - _tmp54 * _tmp76);
  const Scalar _tmp126 = Scalar(1.4083112389913199) * _tmp125;
  const Scalar _tmp127 = std::pow(_tmp125, Scalar(-2));
  const Scalar _tmp128 = _tmp37 + _tmp97;
  const Scalar _tmp129 = _tmp100 * _tmp128;
  const Scalar _tmp130 = Scalar(1.0) / (-_tmp129 - _tmp41 + _tmp99);
  const Scalar _tmp131 = Scalar(1.0) * _tmp130;
  const Scalar _tmp132 = _tmp128 * _tmp130;
  const Scalar _tmp133 = _tmp114 + _tmp116 * _tmp132 - _tmp117 * _tmp95;
  const Scalar _tmp134 = Scalar(1.0) * _tmp98;
  const Scalar _tmp135 = fh1 * (_tmp88 + _tmp91);
  const Scalar _tmp136 = -_tmp135 * _tmp80 - Scalar(40.024799999999999) * _tmp18 - _tmp22 * fv1;
  const Scalar _tmp137 = _tmp98 * (_tmp129 * _tmp131 + Scalar(1.0));
  const Scalar _tmp138 = _tmp100 * _tmp131;
  const Scalar _tmp139 = -Scalar(1.0) * _tmp137 + Scalar(1.0) * _tmp138;
  const Scalar _tmp140 = -_tmp104 + _tmp108 * _tmp132 - _tmp109 * _tmp95;
  const Scalar _tmp141 = _tmp135 * _tmp81 + Scalar(40.024799999999999) * _tmp29 + _tmp33 * fv1;
  const Scalar _tmp142 = _tmp128 * _tmp98;
  const Scalar _tmp143 = Scalar(1.0) * _tmp131 * _tmp142 - Scalar(1.0) * _tmp131;
  const Scalar _tmp144 = _tmp121 * _tmp130;
  const Scalar _tmp145 = _tmp98 * (-_tmp122 * _tmp95 - _tmp128 * _tmp144 + _tmp97);
  const Scalar _tmp146 = -Scalar(1.0) * _tmp144 - Scalar(1.0) * _tmp145 + Scalar(1.0);
  const Scalar _tmp147 = _tmp107 * _tmp131;
  const Scalar _tmp148 = _tmp128 * _tmp147 - _tmp74 * _tmp95;
  const Scalar _tmp149 = Scalar(1.0) * _tmp112 * (_tmp108 * _tmp131 - _tmp134 * _tmp140) +
                         Scalar(1.0) * _tmp119 * (_tmp116 * _tmp131 - _tmp133 * _tmp134) +
                         _tmp120 * _tmp146 + _tmp136 * _tmp139 + _tmp141 * _tmp143 +
                         Scalar(1.0) * _tmp82 * (-_tmp134 * _tmp148 + _tmp147);
  const Scalar _tmp150 = Scalar(1.0) / (_tmp125);
  const Scalar _tmp151 = _tmp19 + _tmp36;
  const Scalar _tmp152 = _tmp124 * _tmp127;
  const Scalar _tmp153 =
      (-_tmp149 * _tmp152 + _tmp150 * (_tmp139 * _tmp151 + _tmp143 * _tmp33 - _tmp146)) /
      std::sqrt(Scalar(_tmp127 * std::pow(_tmp149, Scalar(2)) + 1));
  const Scalar _tmp154 = std::asinh(_tmp149 * _tmp150);
  const Scalar _tmp155 = Scalar(1.4083112389913199) * _tmp124;
  const Scalar _tmp156 = Scalar(0.71007031138673404) * _tmp150;
  const Scalar _tmp157 =
      -_tmp126 * _tmp154 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(Scalar(0.77439999999999998) *
                         std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp56), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp59 - 1), Scalar(2))));
  const Scalar _tmp158 = _tmp156 * _tmp157;
  const Scalar _tmp159 = Scalar(1.0) * _tmp154;
  const Scalar _tmp160 = _tmp131 * _tmp141;
  const Scalar _tmp161 = _tmp112 * _tmp140 * _tmp98 + _tmp119 * _tmp133 * _tmp98 +
                         _tmp120 * _tmp145 + _tmp136 * _tmp137 - _tmp142 * _tmp160 +
                         _tmp148 * _tmp82 * _tmp98;
  const Scalar _tmp162 = _tmp120 * _tmp122;
  const Scalar _tmp163 =
      _tmp110 * _tmp112 * _tmp69 + _tmp118 * _tmp119 * _tmp69 - _tmp162 * _tmp71 - _tmp76 * _tmp82;
  const Scalar _tmp164 = Scalar(1.0) / (_tmp163);
  const Scalar _tmp165 = std::asinh(_tmp161 * _tmp164);
  const Scalar _tmp166 = Scalar(1.4083112389913199) * _tmp122;
  const Scalar _tmp167 = _tmp166 * _tmp71;
  const Scalar _tmp168 = Scalar(1.4083112389913199) * _tmp163;
  const Scalar _tmp169 = _tmp131 * _tmp33;
  const Scalar _tmp170 = std::pow(_tmp163, Scalar(-2));
  const Scalar _tmp171 = _tmp123 * _tmp170;
  const Scalar _tmp172 =
      (-_tmp161 * _tmp171 + _tmp164 * (_tmp137 * _tmp151 - _tmp142 * _tmp169 - _tmp145)) /
      std::sqrt(Scalar(std::pow(_tmp161, Scalar(2)) * _tmp170 + 1));
  const Scalar _tmp173 = Scalar(0.71007031138673404) * _tmp164;
  const Scalar _tmp174 =
      -_tmp165 * _tmp168 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp51), Scalar(2)) +
                     Scalar(0.77439999999999998) *
                         std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp48 - 1), Scalar(2))));
  const Scalar _tmp175 = _tmp173 * _tmp174;
  const Scalar _tmp176 = Scalar(1.0) * _tmp165;
  const Scalar _tmp177 = _tmp109 * _tmp112 + _tmp117 * _tmp119 + _tmp162 + _tmp74 * _tmp82;
  const Scalar _tmp178 = Scalar(1.0) / (_tmp177);
  const Scalar _tmp179 = -_tmp108 * _tmp112 * _tmp130 - _tmp116 * _tmp119 * _tmp130 +
                         _tmp120 * _tmp144 - _tmp136 * _tmp138 - _tmp147 * _tmp82 + _tmp160;
  const Scalar _tmp180 = std::asinh(_tmp178 * _tmp179);
  const Scalar _tmp181 = Scalar(1.0) * _tmp180;
  const Scalar _tmp182 = Scalar(1.4083112389913199) * _tmp177;
  const Scalar _tmp183 =
      -_tmp180 * _tmp182 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(Scalar(0.77439999999999998) *
                         std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp38), Scalar(2)) +
                     std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp42), Scalar(2))));
  const Scalar _tmp184 = Scalar(0.71007031138673404) * _tmp178;
  const Scalar _tmp185 = _tmp183 * _tmp184;
  const Scalar _tmp186 = std::pow(_tmp177, Scalar(-2));
  const Scalar _tmp187 = _tmp122 * _tmp186;
  const Scalar _tmp188 = (_tmp178 * (-_tmp138 * _tmp151 - _tmp144 + _tmp169) + _tmp179 * _tmp187) /
                         std::sqrt(Scalar(std::pow(_tmp179, Scalar(2)) * _tmp186 + 1));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      _tmp35 *
      (-_tmp2 * std::cosh(Scalar(1.0) * _tmp1) +
       _tmp2 * std::cosh(Scalar(0.71007031138673404) * _tmp0 *
                         (-_tmp1 * _tmp35 -
                          Scalar(125.0) *
                              std::sqrt(Scalar(
                                  Scalar(0.77439999999999998) *
                                      std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp23 - 1),
                                               Scalar(2)) +
                                  std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp34 - 1),
                                           Scalar(2)))))));
  _res(1, 0) = _tmp126 * (-Scalar(1.0) * _tmp153 * std::cosh(_tmp159) -
                          (-Scalar(0.71007031138673404) * _tmp152 * _tmp157 +
                           _tmp156 * (-_tmp126 * _tmp153 - _tmp154 * _tmp155)) *
                              std::cosh(_tmp158)) +
               _tmp155 * (-std::sinh(_tmp158) - std::sinh(_tmp159));
  _res(2, 0) = _tmp167 * (-std::sinh(_tmp175) - std::sinh(_tmp176)) +
               _tmp168 * (-Scalar(1.0) * _tmp172 * std::cosh(_tmp176) -
                          (-Scalar(0.71007031138673404) * _tmp171 * _tmp174 +
                           _tmp173 * (-_tmp165 * _tmp167 - _tmp168 * _tmp172)) *
                              std::cosh(_tmp175));
  _res(3, 0) = -_tmp166 * (-std::sinh(_tmp181) - std::sinh(_tmp185)) +
               _tmp182 * (-Scalar(1.0) * _tmp188 * std::cosh(_tmp181) -
                          (Scalar(0.71007031138673404) * _tmp183 * _tmp187 +
                           _tmp184 * (_tmp166 * _tmp180 - _tmp182 * _tmp188)) *
                              std::cosh(_tmp185));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
