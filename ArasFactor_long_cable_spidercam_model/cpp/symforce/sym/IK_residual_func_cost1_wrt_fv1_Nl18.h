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
 * Symbolic function: IK_residual_func_cost1_wrt_fv1_Nl18
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1WrtFv1Nl18(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 609

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (193)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = std::asinh(_tmp0 * fv1);
  const Scalar _tmp2 = Scalar(1.0) * _tmp0 /
                       std::sqrt(Scalar(1 + std::pow(fv1, Scalar(2)) / std::pow(fh1, Scalar(2))));
  const Scalar _tmp3 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp4 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp6 = 1 - 2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp7 = Scalar(0.20999999999999999) * _tmp4 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp8 = -_tmp7;
  const Scalar _tmp9 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp10 = 2 * _tmp5 * _tmp9;
  const Scalar _tmp11 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                        2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp12 = _tmp11 * _tmp3;
  const Scalar _tmp13 = _tmp10 + _tmp12;
  const Scalar _tmp14 = -Scalar(0.010999999999999999) * _tmp13;
  const Scalar _tmp15 = 2 * _tmp3;
  const Scalar _tmp16 = _tmp15 * _tmp9;
  const Scalar _tmp17 = _tmp11 * _tmp5;
  const Scalar _tmp18 = Scalar(0.20999999999999999) * _tmp16 - Scalar(0.20999999999999999) * _tmp17;
  const Scalar _tmp19 = _tmp14 + _tmp18;
  const Scalar _tmp20 = _tmp19 + _tmp8;
  const Scalar _tmp21 = _tmp20 + position_vector(0, 0);
  const Scalar _tmp22 = Scalar(0.20999999999999999) * _tmp16 + Scalar(0.20999999999999999) * _tmp17;
  const Scalar _tmp23 = -_tmp22;
  const Scalar _tmp24 = _tmp15 * _tmp5;
  const Scalar _tmp25 = _tmp11 * _tmp9;
  const Scalar _tmp26 = _tmp24 - _tmp25;
  const Scalar _tmp27 = Scalar(0.010999999999999999) * _tmp26;
  const Scalar _tmp28 = -_tmp27;
  const Scalar _tmp29 = -2 * std::pow(_tmp9, Scalar(2));
  const Scalar _tmp30 = Scalar(0.20999999999999999) * _tmp29 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp31 = _tmp28 + _tmp30;
  const Scalar _tmp32 = _tmp23 + _tmp31;
  const Scalar _tmp33 = _tmp32 + position_vector(1, 0);
  const Scalar _tmp34 = Scalar(1.4083112389913199) * fh1;
  const Scalar _tmp35 = -_tmp30;
  const Scalar _tmp36 = _tmp28 + _tmp35;
  const Scalar _tmp37 = _tmp23 + _tmp36;
  const Scalar _tmp38 = _tmp37 + position_vector(1, 0);
  const Scalar _tmp39 = _tmp38 + Scalar(110.0);
  const Scalar _tmp40 = _tmp14 - _tmp18;
  const Scalar _tmp41 = _tmp40 + _tmp8;
  const Scalar _tmp42 = _tmp41 + position_vector(0, 0);
  const Scalar _tmp43 = _tmp42 + Scalar(125.0);
  const Scalar _tmp44 = Scalar(1.0) / (_tmp43);
  const Scalar _tmp45 = _tmp39 * _tmp44;
  const Scalar _tmp46 = _tmp22 + _tmp36;
  const Scalar _tmp47 = _tmp46 + position_vector(1, 0);
  const Scalar _tmp48 = _tmp47 + Scalar(110.0);
  const Scalar _tmp49 = _tmp40 + _tmp7;
  const Scalar _tmp50 = _tmp49 + position_vector(0, 0);
  const Scalar _tmp51 = _tmp50 + Scalar(-125.0);
  const Scalar _tmp52 = std::pow(Scalar(std::pow(_tmp48, Scalar(2)) + std::pow(_tmp51, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp53 = _tmp51 * _tmp52;
  const Scalar _tmp54 =
      std::sqrt(Scalar(std::pow(_tmp39, Scalar(2)) + std::pow(_tmp43, Scalar(2))));
  const Scalar _tmp55 = Scalar(1.0) / (_tmp54);
  const Scalar _tmp56 = _tmp44 * _tmp54;
  const Scalar _tmp57 = _tmp56 * (-_tmp37 * _tmp43 * _tmp55 + _tmp39 * _tmp41 * _tmp55);
  const Scalar _tmp58 = _tmp48 * _tmp52;
  const Scalar _tmp59 = _tmp46 * _tmp53 - _tmp49 * _tmp58 + _tmp53 * _tmp57;
  const Scalar _tmp60 = Scalar(1.0) / (_tmp45 * _tmp53 - _tmp58);
  const Scalar _tmp61 = _tmp45 * _tmp60;
  const Scalar _tmp62 = Scalar(0.20999999999999999) * _tmp10 - Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp63 = -Scalar(0.010999999999999999) * _tmp29 -
                        Scalar(0.010999999999999999) * _tmp4 + Scalar(-0.010999999999999999);
  const Scalar _tmp64 = Scalar(0.20999999999999999) * _tmp24 + Scalar(0.20999999999999999) * _tmp25;
  const Scalar _tmp65 = _tmp63 - _tmp64;
  const Scalar _tmp66 = _tmp62 + _tmp65;
  const Scalar _tmp67 = -_tmp62;
  const Scalar _tmp68 = _tmp65 + _tmp67;
  const Scalar _tmp69 = _tmp53 * _tmp68;
  const Scalar _tmp70 = -_tmp53 * _tmp66 + _tmp69;
  const Scalar _tmp71 = _tmp60 * (-_tmp45 * _tmp69 + _tmp58 * _tmp66);
  const Scalar _tmp72 = _tmp45 * _tmp68 + _tmp45 * _tmp71;
  const Scalar _tmp73 = Scalar(1.0) * _tmp37;
  const Scalar _tmp74 = -_tmp73;
  const Scalar _tmp75 = Scalar(1.0) / (_tmp46 + _tmp74);
  const Scalar _tmp76 = Scalar(1.0) * _tmp41;
  const Scalar _tmp77 = _tmp75 * (-_tmp49 + _tmp76);
  const Scalar _tmp78 = _tmp61 * _tmp70 - _tmp68 - _tmp72 * _tmp77;
  const Scalar _tmp79 = _tmp22 + _tmp31;
  const Scalar _tmp80 = _tmp79 + position_vector(1, 0);
  const Scalar _tmp81 = _tmp80 + Scalar(-110.0);
  const Scalar _tmp82 = _tmp19 + _tmp7;
  const Scalar _tmp83 = _tmp82 + position_vector(0, 0);
  const Scalar _tmp84 = _tmp83 + Scalar(-125.0);
  const Scalar _tmp85 = std::pow(Scalar(std::pow(_tmp81, Scalar(2)) + std::pow(_tmp84, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp86 = _tmp84 * _tmp85;
  const Scalar _tmp87 = _tmp68 * _tmp86;
  const Scalar _tmp88 = _tmp63 + _tmp64;
  const Scalar _tmp89 = _tmp62 + _tmp88;
  const Scalar _tmp90 = _tmp81 * _tmp85;
  const Scalar _tmp91 = _tmp45 * _tmp86 - _tmp90;
  const Scalar _tmp92 = -_tmp45 * _tmp87 - _tmp71 * _tmp91 + _tmp89 * _tmp90;
  const Scalar _tmp93 = _tmp60 * _tmp91;
  const Scalar _tmp94 = -_tmp70 * _tmp93 - _tmp77 * _tmp92 - _tmp86 * _tmp89 + _tmp87;
  const Scalar _tmp95 = Scalar(1.0) / (_tmp94);
  const Scalar _tmp96 = _tmp57 * _tmp86 - _tmp59 * _tmp93 + _tmp79 * _tmp86 - _tmp82 * _tmp90;
  const Scalar _tmp97 = _tmp95 * _tmp96;
  const Scalar _tmp98 = Scalar(1.0) / (_tmp96);
  const Scalar _tmp99 = _tmp94 * _tmp98;
  const Scalar _tmp100 = _tmp99 * (-_tmp57 + _tmp59 * _tmp61 - _tmp78 * _tmp97);
  const Scalar _tmp101 = _tmp100 + _tmp78;
  const Scalar _tmp102 = _tmp91 * _tmp95;
  const Scalar _tmp103 = -_tmp101 * _tmp102 - _tmp45;
  const Scalar _tmp104 = _tmp53 * _tmp60;
  const Scalar _tmp105 = _tmp86 * _tmp95;
  const Scalar _tmp106 = _tmp21 + Scalar(125.0);
  const Scalar _tmp107 = _tmp33 + Scalar(-110.0);
  const Scalar _tmp108 =
      std::pow(Scalar(std::pow(_tmp106, Scalar(2)) + std::pow(_tmp107, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp109 = _tmp106 * _tmp108;
  const Scalar _tmp110 = _tmp109 * fh1;
  const Scalar _tmp111 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp112 = _tmp73 * _tmp77 + _tmp76;
  const Scalar _tmp113 = 0;
  const Scalar _tmp114 = _tmp113 * _tmp93;
  const Scalar _tmp115 = _tmp56 * (_tmp113 * _tmp86 - _tmp114 * _tmp53);
  const Scalar _tmp116 = Scalar(1.0) * _tmp98;
  const Scalar _tmp117 = Scalar(1.0) * _tmp60;
  const Scalar _tmp118 = _tmp117 * _tmp91 * _tmp98;
  const Scalar _tmp119 = _tmp107 * _tmp108;
  const Scalar _tmp120 = fh1 * (-_tmp109 * _tmp32 + _tmp119 * _tmp20);
  const Scalar _tmp121 = Scalar(1.0) * _tmp71;
  const Scalar _tmp122 = -_tmp117 * _tmp70 + _tmp121 * _tmp77;
  const Scalar _tmp123 = _tmp99 * (-_tmp117 * _tmp59 - _tmp122 * _tmp97);
  const Scalar _tmp124 = _tmp122 + _tmp123;
  const Scalar _tmp125 = -_tmp102 * _tmp124 + Scalar(1.0);
  const Scalar _tmp126 = _tmp119 * fh1;
  const Scalar _tmp127 = -_tmp110 * _tmp56 * (_tmp101 * _tmp105 + _tmp103 * _tmp104 + Scalar(1.0)) -
                         _tmp111 * _tmp115 -
                         _tmp120 * _tmp56 * (_tmp116 * _tmp86 - _tmp118 * _tmp53) -
                         _tmp126 * _tmp56 * (_tmp104 * _tmp125 + _tmp105 * _tmp124);
  const Scalar _tmp128 = std::pow(_tmp127, Scalar(-2));
  const Scalar _tmp129 = _tmp115 * _tmp128;
  const Scalar _tmp130 = _tmp74 + _tmp79;
  const Scalar _tmp131 = _tmp130 * _tmp77;
  const Scalar _tmp132 = Scalar(1.0) / (-_tmp131 + _tmp76 - _tmp82);
  const Scalar _tmp133 = _tmp112 * _tmp132;
  const Scalar _tmp134 = _tmp75 * (-_tmp113 * _tmp92 - _tmp130 * _tmp133 + _tmp74);
  const Scalar _tmp135 = Scalar(1.0) * _tmp132;
  const Scalar _tmp136 = -Scalar(1.0) * _tmp112 * _tmp135 - Scalar(1.0) * _tmp134 + Scalar(1.0);
  const Scalar _tmp137 = fh1 * (_tmp67 + _tmp88);
  const Scalar _tmp138 = -_tmp119 * _tmp137 - Scalar(40.024799999999999) * _tmp26 - _tmp32 * fv1;
  const Scalar _tmp139 = _tmp75 * (_tmp131 * _tmp135 + Scalar(1.0));
  const Scalar _tmp140 = _tmp135 * _tmp77;
  const Scalar _tmp141 = -Scalar(1.0) * _tmp139 + Scalar(1.0) * _tmp140;
  const Scalar _tmp142 = _tmp109 * _tmp137 + Scalar(40.024799999999999) * _tmp13 + _tmp20 * fv1;
  const Scalar _tmp143 = _tmp130 * _tmp75;
  const Scalar _tmp144 = Scalar(1.0) * _tmp135 * _tmp143 - Scalar(1.0) * _tmp135;
  const Scalar _tmp145 = _tmp92 * _tmp95;
  const Scalar _tmp146 = _tmp130 * _tmp132;
  const Scalar _tmp147 = _tmp100 * _tmp146 - _tmp101 * _tmp145 + _tmp72;
  const Scalar _tmp148 = Scalar(1.0) * _tmp75;
  const Scalar _tmp149 = Scalar(1.0) * fh1;
  const Scalar _tmp150 = -_tmp121 + _tmp123 * _tmp146 - _tmp124 * _tmp145;
  const Scalar _tmp151 = _tmp135 * _tmp99;
  const Scalar _tmp152 = -_tmp116 * _tmp92 + _tmp130 * _tmp151;
  const Scalar _tmp153 = _tmp109 * _tmp149 * (_tmp100 * _tmp135 - _tmp147 * _tmp148) +
                         _tmp111 * _tmp136 +
                         _tmp119 * _tmp149 * (_tmp123 * _tmp135 - _tmp148 * _tmp150) +
                         Scalar(1.0) * _tmp120 * (-_tmp148 * _tmp152 + _tmp151) +
                         _tmp138 * _tmp141 + _tmp142 * _tmp144;
  const Scalar _tmp154 = Scalar(1.0) / (_tmp127);
  const Scalar _tmp155 = std::asinh(_tmp153 * _tmp154);
  const Scalar _tmp156 = Scalar(1.0) * _tmp155;
  const Scalar _tmp157 = _tmp22 + _tmp27 + _tmp35;
  const Scalar _tmp158 =
      (-_tmp129 * _tmp153 + _tmp154 * (-_tmp136 + _tmp141 * _tmp157 + _tmp144 * _tmp20)) /
      std::sqrt(Scalar(_tmp128 * std::pow(_tmp153, Scalar(2)) + 1));
  const Scalar _tmp159 = Scalar(1.4083112389913199) * _tmp127;
  const Scalar _tmp160 =
      -_tmp155 * _tmp159 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(Scalar(0.77439999999999998) *
                         std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp38 - 1), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp42 - 1), Scalar(2))));
  const Scalar _tmp161 = Scalar(1.4083112389913199) * _tmp115;
  const Scalar _tmp162 = Scalar(0.71007031138673404) * _tmp154;
  const Scalar _tmp163 = _tmp160 * _tmp162;
  const Scalar _tmp164 = _tmp135 * _tmp142;
  const Scalar _tmp165 = _tmp110 * _tmp147 * _tmp75 + _tmp111 * _tmp134 +
                         _tmp120 * _tmp152 * _tmp75 + _tmp126 * _tmp150 * _tmp75 +
                         _tmp138 * _tmp139 - _tmp143 * _tmp164;
  const Scalar _tmp166 = _tmp111 * _tmp113;
  const Scalar _tmp167 = _tmp103 * _tmp110 * _tmp60 - _tmp118 * _tmp120 +
                         _tmp125 * _tmp126 * _tmp60 - _tmp166 * _tmp93;
  const Scalar _tmp168 = Scalar(1.0) / (_tmp167);
  const Scalar _tmp169 = std::asinh(_tmp165 * _tmp168);
  const Scalar _tmp170 = Scalar(1.0) * _tmp169;
  const Scalar _tmp171 = Scalar(1.4083112389913199) * _tmp167;
  const Scalar _tmp172 =
      -_tmp169 * _tmp171 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp50), Scalar(2)) +
                     Scalar(0.77439999999999998) *
                         std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp47 - 1), Scalar(2))));
  const Scalar _tmp173 = Scalar(0.71007031138673404) * _tmp168;
  const Scalar _tmp174 = _tmp172 * _tmp173;
  const Scalar _tmp175 = Scalar(1.4083112389913199) * _tmp113;
  const Scalar _tmp176 = _tmp175 * _tmp93;
  const Scalar _tmp177 = std::pow(_tmp167, Scalar(-2));
  const Scalar _tmp178 = _tmp114 * _tmp177;
  const Scalar _tmp179 = _tmp135 * _tmp20;
  const Scalar _tmp180 =
      (-_tmp165 * _tmp178 + _tmp168 * (-_tmp134 + _tmp139 * _tmp157 - _tmp143 * _tmp179)) /
      std::sqrt(Scalar(std::pow(_tmp165, Scalar(2)) * _tmp177 + 1));
  const Scalar _tmp181 =
      _tmp101 * _tmp110 * _tmp95 + _tmp116 * _tmp120 + _tmp124 * _tmp126 * _tmp95 + _tmp166;
  const Scalar _tmp182 = std::pow(_tmp181, Scalar(-2));
  const Scalar _tmp183 = _tmp113 * _tmp182;
  const Scalar _tmp184 = Scalar(1.0) / (_tmp181);
  const Scalar _tmp185 = -_tmp100 * _tmp110 * _tmp132 + _tmp111 * _tmp133 - _tmp120 * _tmp151 -
                         _tmp123 * _tmp126 * _tmp132 - _tmp138 * _tmp140 + _tmp164;
  const Scalar _tmp186 = std::asinh(_tmp184 * _tmp185);
  const Scalar _tmp187 = Scalar(1.4083112389913199) * _tmp181;
  const Scalar _tmp188 =
      -_tmp186 * _tmp187 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(Scalar(0.77439999999999998) *
                         std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp80), Scalar(2)) +
                     std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp83), Scalar(2))));
  const Scalar _tmp189 = (_tmp183 * _tmp185 + _tmp184 * (-_tmp133 - _tmp140 * _tmp157 + _tmp179)) /
                         std::sqrt(Scalar(_tmp182 * std::pow(_tmp185, Scalar(2)) + 1));
  const Scalar _tmp190 = Scalar(0.71007031138673404) * _tmp184;
  const Scalar _tmp191 = _tmp188 * _tmp190;
  const Scalar _tmp192 = Scalar(1.0) * _tmp186;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      -_tmp34 *
      (_tmp2 * std::sinh(Scalar(1.0) * _tmp1) +
       _tmp2 * std::sinh(Scalar(0.71007031138673404) * _tmp0 *
                         (-_tmp1 * _tmp34 -
                          Scalar(125.0) *
                              std::sqrt(Scalar(
                                  Scalar(0.77439999999999998) *
                                      std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp33),
                                               Scalar(2)) +
                                  std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp21 - 1),
                                           Scalar(2)))))));
  _res(1, 0) =
      -_tmp159 *
          (-Scalar(34.083374946563197) * _tmp129 + Scalar(1.0) * _tmp158 * std::sinh(_tmp156) -
           (-Scalar(0.71007031138673404) * _tmp129 * _tmp160 +
            _tmp162 * (-_tmp155 * _tmp161 - _tmp158 * _tmp159)) *
               std::sinh(_tmp163)) -
      _tmp161 * (Scalar(34.083374946563197) * _tmp154 + std::cosh(_tmp156) - std::cosh(_tmp163));
  _res(2, 0) =
      -_tmp171 *
          (-Scalar(34.083374946563197) * _tmp178 + Scalar(1.0) * _tmp180 * std::sinh(_tmp170) -
           (-Scalar(0.71007031138673404) * _tmp172 * _tmp178 +
            _tmp173 * (-_tmp169 * _tmp176 - _tmp171 * _tmp180)) *
               std::sinh(_tmp174)) -
      _tmp176 * (Scalar(34.083374946563197) * _tmp168 + std::cosh(_tmp170) - std::cosh(_tmp174));
  _res(3, 0) =
      _tmp175 * (Scalar(34.083374946563197) * _tmp184 - std::cosh(_tmp191) + std::cosh(_tmp192)) -
      _tmp187 * (Scalar(34.083374946563197) * _tmp183 + Scalar(1.0) * _tmp189 * std::sinh(_tmp192) -
                 (Scalar(0.71007031138673404) * _tmp183 * _tmp188 +
                  _tmp190 * (_tmp175 * _tmp186 - _tmp187 * _tmp189)) *
                     std::sinh(_tmp191));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
