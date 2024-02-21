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
 * Symbolic function: IK_residual_func_cost1_wrt_fh1_Nl3
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1WrtFh1Nl3(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 660

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (214)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp2 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp3 = 2 * _tmp1 * _tmp2;
  const Scalar _tmp4 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp5 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp6 = _tmp4 * _tmp5;
  const Scalar _tmp7 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp8 = -_tmp7;
  const Scalar _tmp9 = 2 * _tmp4;
  const Scalar _tmp10 = _tmp1 * _tmp9;
  const Scalar _tmp11 = _tmp2 * _tmp5;
  const Scalar _tmp12 = _tmp10 - _tmp11;
  const Scalar _tmp13 = -Scalar(0.010999999999999999) * _tmp12;
  const Scalar _tmp14 = -2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp15 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp16 = Scalar(0.20999999999999999) * _tmp14 +
                        Scalar(0.20999999999999999) * _tmp15 + Scalar(0.20999999999999999);
  const Scalar _tmp17 = _tmp13 - _tmp16;
  const Scalar _tmp18 = _tmp17 + _tmp8;
  const Scalar _tmp19 = _tmp18 + position_vector(1, 0);
  const Scalar _tmp20 = 1 - 2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp15 + Scalar(0.20999999999999999) * _tmp20;
  const Scalar _tmp22 = -_tmp21;
  const Scalar _tmp23 = _tmp2 * _tmp9;
  const Scalar _tmp24 = _tmp1 * _tmp5;
  const Scalar _tmp25 = _tmp23 + _tmp24;
  const Scalar _tmp26 = -Scalar(0.010999999999999999) * _tmp25;
  const Scalar _tmp27 = Scalar(0.20999999999999999) * _tmp3 - Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp28 = _tmp26 - _tmp27;
  const Scalar _tmp29 = _tmp22 + _tmp28;
  const Scalar _tmp30 = _tmp29 + position_vector(0, 0);
  const Scalar _tmp31 = _tmp0 * fv1;
  const Scalar _tmp32 = std::asinh(_tmp31);
  const Scalar _tmp33 = Scalar(1.4083112389913199) * _tmp32;
  const Scalar _tmp34 =
      -_tmp33 * fh1 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(Scalar(0.77439999999999998) *
                         std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp19 - 1), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp30 - 1), Scalar(2))));
  const Scalar _tmp35 = Scalar(0.71007031138673404) * _tmp0;
  const Scalar _tmp36 = _tmp34 * _tmp35;
  const Scalar _tmp37 = Scalar(1.0) * _tmp32;
  const Scalar _tmp38 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp39 =
      std::pow(Scalar(_tmp38 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp40 = _tmp26 + _tmp27;
  const Scalar _tmp41 = _tmp21 + _tmp40;
  const Scalar _tmp42 = _tmp41 + position_vector(0, 0);
  const Scalar _tmp43 = _tmp13 + _tmp16;
  const Scalar _tmp44 = _tmp43 + _tmp7;
  const Scalar _tmp45 = _tmp44 + position_vector(1, 0);
  const Scalar _tmp46 = _tmp45 + Scalar(-110.0);
  const Scalar _tmp47 = _tmp42 + Scalar(-125.0);
  const Scalar _tmp48 =
      std::sqrt(Scalar(std::pow(_tmp46, Scalar(2)) + std::pow(_tmp47, Scalar(2))));
  const Scalar _tmp49 = Scalar(1.0) / (_tmp48);
  const Scalar _tmp50 = Scalar(1.0) / (_tmp47);
  const Scalar _tmp51 = _tmp48 * _tmp50;
  const Scalar _tmp52 = _tmp51 * (_tmp41 * _tmp46 * _tmp49 - _tmp44 * _tmp47 * _tmp49);
  const Scalar _tmp53 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp11;
  const Scalar _tmp54 =
      -Scalar(0.010999999999999999) * _tmp14 - Scalar(0.010999999999999999) * _tmp20;
  const Scalar _tmp55 = Scalar(0.20999999999999999) * _tmp23 - Scalar(0.20999999999999999) * _tmp24;
  const Scalar _tmp56 = _tmp54 - _tmp55;
  const Scalar _tmp57 = _tmp53 + _tmp56;
  const Scalar _tmp58 = _tmp22 + _tmp40;
  const Scalar _tmp59 = _tmp58 + position_vector(0, 0);
  const Scalar _tmp60 = _tmp59 + Scalar(125.0);
  const Scalar _tmp61 = _tmp43 + _tmp8;
  const Scalar _tmp62 = _tmp61 + position_vector(1, 0);
  const Scalar _tmp63 = _tmp62 + Scalar(-110.0);
  const Scalar _tmp64 = std::pow(Scalar(std::pow(_tmp60, Scalar(2)) + std::pow(_tmp63, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp65 = _tmp60 * _tmp64;
  const Scalar _tmp66 = _tmp54 + _tmp55;
  const Scalar _tmp67 = _tmp53 + _tmp66;
  const Scalar _tmp68 = -_tmp57 * _tmp65 + _tmp65 * _tmp67;
  const Scalar _tmp69 = _tmp63 * _tmp64;
  const Scalar _tmp70 = _tmp46 * _tmp50;
  const Scalar _tmp71 = Scalar(1.0) / (_tmp65 * _tmp70 - _tmp69);
  const Scalar _tmp72 = _tmp70 * _tmp71;
  const Scalar _tmp73 = _tmp67 * _tmp70;
  const Scalar _tmp74 = _tmp57 * _tmp69 - _tmp65 * _tmp73;
  const Scalar _tmp75 = _tmp72 * _tmp74 + _tmp73;
  const Scalar _tmp76 = Scalar(1.0) * _tmp44;
  const Scalar _tmp77 = -_tmp76;
  const Scalar _tmp78 = Scalar(1.0) / (_tmp61 + _tmp77);
  const Scalar _tmp79 = Scalar(1.0) * _tmp41;
  const Scalar _tmp80 = -_tmp58 + _tmp79;
  const Scalar _tmp81 = _tmp78 * _tmp80;
  const Scalar _tmp82 = -_tmp67 + _tmp68 * _tmp72 - _tmp75 * _tmp81;
  const Scalar _tmp83 = _tmp17 + _tmp7;
  const Scalar _tmp84 = _tmp83 + position_vector(1, 0);
  const Scalar _tmp85 = _tmp84 + Scalar(110.0);
  const Scalar _tmp86 = _tmp21 + _tmp28;
  const Scalar _tmp87 = _tmp86 + position_vector(0, 0);
  const Scalar _tmp88 = _tmp87 + Scalar(-125.0);
  const Scalar _tmp89 = std::pow(Scalar(std::pow(_tmp85, Scalar(2)) + std::pow(_tmp88, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp90 = _tmp85 * _tmp89;
  const Scalar _tmp91 = _tmp88 * _tmp89;
  const Scalar _tmp92 = _tmp70 * _tmp91 - _tmp90;
  const Scalar _tmp93 = _tmp71 * _tmp92;
  const Scalar _tmp94 = -_tmp53;
  const Scalar _tmp95 = _tmp66 + _tmp94;
  const Scalar _tmp96 = -_tmp73 * _tmp91 - _tmp74 * _tmp93 + _tmp90 * _tmp95;
  const Scalar _tmp97 = _tmp67 * _tmp91 - _tmp68 * _tmp93 - _tmp81 * _tmp96 - _tmp91 * _tmp95;
  const Scalar _tmp98 = Scalar(1.0) / (_tmp97);
  const Scalar _tmp99 = _tmp52 * _tmp65 - _tmp58 * _tmp69 + _tmp61 * _tmp65;
  const Scalar _tmp100 = _tmp52 * _tmp91 + _tmp83 * _tmp91 - _tmp86 * _tmp90 - _tmp93 * _tmp99;
  const Scalar _tmp101 = _tmp100 * _tmp98;
  const Scalar _tmp102 = Scalar(1.0) / (_tmp100);
  const Scalar _tmp103 = _tmp102 * _tmp97;
  const Scalar _tmp104 = _tmp103 * (-_tmp101 * _tmp82 - _tmp52 + _tmp72 * _tmp99);
  const Scalar _tmp105 = _tmp104 + _tmp82;
  const Scalar _tmp106 = _tmp91 * _tmp98;
  const Scalar _tmp107 = _tmp92 * _tmp98;
  const Scalar _tmp108 = _tmp71 * (-_tmp105 * _tmp107 - _tmp70);
  const Scalar _tmp109 = _tmp30 + Scalar(125.0);
  const Scalar _tmp110 = _tmp19 + Scalar(110.0);
  const Scalar _tmp111 =
      std::pow(Scalar(std::pow(_tmp109, Scalar(2)) + std::pow(_tmp110, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp112 = _tmp109 * _tmp111;
  const Scalar _tmp113 = _tmp112 * _tmp51 * (_tmp105 * _tmp106 + _tmp108 * _tmp65 + Scalar(1.0));
  const Scalar _tmp114 = _tmp110 * _tmp111;
  const Scalar _tmp115 = -_tmp112 * _tmp18 + _tmp114 * _tmp29;
  const Scalar _tmp116 = Scalar(1.0) * _tmp102;
  const Scalar _tmp117 = Scalar(1.0) * _tmp71;
  const Scalar _tmp118 = _tmp102 * _tmp117 * _tmp92;
  const Scalar _tmp119 = _tmp115 * _tmp51 * (_tmp116 * _tmp91 - _tmp118 * _tmp65);
  const Scalar _tmp120 = Scalar(1.0) * _tmp78;
  const Scalar _tmp121 = -_tmp117 * _tmp68 + _tmp120 * _tmp71 * _tmp74 * _tmp80;
  const Scalar _tmp122 = _tmp103 * (-_tmp101 * _tmp121 - _tmp117 * _tmp99);
  const Scalar _tmp123 = _tmp121 + _tmp122;
  const Scalar _tmp124 = _tmp71 * (-_tmp107 * _tmp123 + Scalar(1.0));
  const Scalar _tmp125 = _tmp114 * _tmp51 * (_tmp106 * _tmp123 + _tmp124 * _tmp65);
  const Scalar _tmp126 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp127 = _tmp76 * _tmp81 + _tmp79;
  const Scalar _tmp128 = 0;
  const Scalar _tmp129 = _tmp128 * _tmp98;
  const Scalar _tmp130 = -_tmp113 * fh1 - _tmp119 * fh1 - _tmp125 * fh1 -
                         _tmp126 * _tmp51 * (-_tmp129 * _tmp65 * _tmp93 + _tmp129 * _tmp91);
  const Scalar _tmp131 = Scalar(1.0) / (_tmp130);
  const Scalar _tmp132 = _tmp77 + _tmp83;
  const Scalar _tmp133 = _tmp132 * _tmp81;
  const Scalar _tmp134 = Scalar(1.0) / (-_tmp133 + _tmp79 - _tmp86);
  const Scalar _tmp135 = Scalar(1.0) * _tmp134;
  const Scalar _tmp136 = _tmp103 * _tmp135;
  const Scalar _tmp137 = -_tmp116 * _tmp96 + _tmp132 * _tmp136;
  const Scalar _tmp138 = Scalar(1.0) * _tmp115 * (-_tmp120 * _tmp137 + _tmp136);
  const Scalar _tmp139 = _tmp96 * _tmp98;
  const Scalar _tmp140 = _tmp132 * _tmp134;
  const Scalar _tmp141 = _tmp104 * _tmp140 - _tmp105 * _tmp139 + _tmp75;
  const Scalar _tmp142 = Scalar(1.0) * _tmp112 * (_tmp104 * _tmp135 - _tmp120 * _tmp141);
  const Scalar _tmp143 = _tmp56 + _tmp94;
  const Scalar _tmp144 = _tmp143 * fh1;
  const Scalar _tmp145 = _tmp112 * _tmp144 + Scalar(40.024799999999999) * _tmp25 + _tmp29 * fv1;
  const Scalar _tmp146 = _tmp132 * _tmp78;
  const Scalar _tmp147 = _tmp135 * _tmp146;
  const Scalar _tmp148 = -Scalar(1.0) * _tmp135 + Scalar(1.0) * _tmp147;
  const Scalar _tmp149 = -_tmp117 * _tmp74 + _tmp122 * _tmp140 - _tmp123 * _tmp139;
  const Scalar _tmp150 = Scalar(1.0) * _tmp114 * (-_tmp120 * _tmp149 + _tmp122 * _tmp135);
  const Scalar _tmp151 = _tmp127 * _tmp134;
  const Scalar _tmp152 = -_tmp128 * _tmp139 - _tmp132 * _tmp151 + _tmp77;
  const Scalar _tmp153 = -_tmp114 * _tmp144 - Scalar(40.024799999999999) * _tmp12 - _tmp18 * fv1;
  const Scalar _tmp154 = _tmp135 * _tmp81;
  const Scalar _tmp155 = _tmp133 * _tmp135 + Scalar(1.0);
  const Scalar _tmp156 = -Scalar(1.0) * _tmp120 * _tmp155 + Scalar(1.0) * _tmp154;
  const Scalar _tmp157 =
      Scalar(1.0) * _tmp126 * (-_tmp120 * _tmp152 - _tmp127 * _tmp135 + Scalar(1.0)) +
      _tmp138 * fh1 + _tmp142 * fh1 + _tmp145 * _tmp148 + _tmp150 * fh1 + _tmp153 * _tmp156;
  const Scalar _tmp158 = std::asinh(_tmp131 * _tmp157);
  const Scalar _tmp159 = Scalar(1.4083112389913199) * _tmp130;
  const Scalar _tmp160 =
      -_tmp158 * _tmp159 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp42), Scalar(2)) +
                     Scalar(0.77439999999999998) *
                         std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp45), Scalar(2))));
  const Scalar _tmp161 = Scalar(0.71007031138673404) * _tmp131;
  const Scalar _tmp162 = _tmp160 * _tmp161;
  const Scalar _tmp163 = -_tmp113 - _tmp119 - _tmp125;
  const Scalar _tmp164 = Scalar(1.4083112389913199) * _tmp163;
  const Scalar _tmp165 = _tmp114 * _tmp143;
  const Scalar _tmp166 = _tmp112 * _tmp143;
  const Scalar _tmp167 = std::pow(_tmp130, Scalar(-2));
  const Scalar _tmp168 = _tmp163 * _tmp167;
  const Scalar _tmp169 =
      (_tmp131 * (_tmp138 + _tmp142 + _tmp148 * _tmp166 + _tmp150 - _tmp156 * _tmp165) -
       _tmp157 * _tmp168) /
      std::sqrt(Scalar(std::pow(_tmp157, Scalar(2)) * _tmp167 + 1));
  const Scalar _tmp170 = Scalar(1.0) * _tmp158;
  const Scalar _tmp171 = _tmp114 * _tmp124;
  const Scalar _tmp172 = _tmp108 * _tmp112;
  const Scalar _tmp173 = _tmp115 * _tmp118;
  const Scalar _tmp174 = _tmp171 + _tmp172 - _tmp173;
  const Scalar _tmp175 = _tmp126 * _tmp129;
  const Scalar _tmp176 = _tmp171 * fh1 + _tmp172 * fh1 - _tmp173 * fh1 - _tmp175 * _tmp93;
  const Scalar _tmp177 = Scalar(1.0) / (_tmp176);
  const Scalar _tmp178 = _tmp115 * _tmp137 * _tmp78;
  const Scalar _tmp179 = _tmp112 * _tmp141 * _tmp78;
  const Scalar _tmp180 = _tmp135 * _tmp145;
  const Scalar _tmp181 = _tmp155 * _tmp78;
  const Scalar _tmp182 = _tmp114 * _tmp149 * _tmp78;
  const Scalar _tmp183 = _tmp126 * _tmp152 * _tmp78 - _tmp146 * _tmp180 + _tmp153 * _tmp181 +
                         _tmp178 * fh1 + _tmp179 * fh1 + _tmp182 * fh1;
  const Scalar _tmp184 = std::asinh(_tmp177 * _tmp183);
  const Scalar _tmp185 = Scalar(1.4083112389913199) * _tmp184;
  const Scalar _tmp186 =
      -_tmp176 * _tmp185 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(Scalar(0.77439999999999998) *
                         std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp62), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp59 - 1), Scalar(2))));
  const Scalar _tmp187 = Scalar(0.71007031138673404) * _tmp177;
  const Scalar _tmp188 = _tmp186 * _tmp187;
  const Scalar _tmp189 = Scalar(1.0) * _tmp184;
  const Scalar _tmp190 = std::pow(_tmp176, Scalar(-2));
  const Scalar _tmp191 = _tmp174 * _tmp190;
  const Scalar _tmp192 =
      (_tmp177 * (-_tmp147 * _tmp166 - _tmp165 * _tmp181 + _tmp178 + _tmp179 + _tmp182) -
       _tmp183 * _tmp191) /
      std::sqrt(Scalar(std::pow(_tmp183, Scalar(2)) * _tmp190 + 1));
  const Scalar _tmp193 = Scalar(1.4083112389913199) * _tmp176;
  const Scalar _tmp194 = _tmp115 * _tmp116;
  const Scalar _tmp195 = _tmp105 * _tmp112 * _tmp98;
  const Scalar _tmp196 = _tmp114 * _tmp123 * _tmp98;
  const Scalar _tmp197 = _tmp175 + _tmp194 * fh1 + _tmp195 * fh1 + _tmp196 * fh1;
  const Scalar _tmp198 = Scalar(1.0) / (_tmp197);
  const Scalar _tmp199 = _tmp115 * _tmp136;
  const Scalar _tmp200 = _tmp114 * _tmp122 * _tmp134;
  const Scalar _tmp201 = _tmp104 * _tmp112 * _tmp134;
  const Scalar _tmp202 = _tmp126 * _tmp151 - _tmp153 * _tmp154 + _tmp180 - _tmp199 * fh1 -
                         _tmp200 * fh1 - _tmp201 * fh1;
  const Scalar _tmp203 = std::asinh(_tmp198 * _tmp202);
  const Scalar _tmp204 = Scalar(1.4083112389913199) * _tmp203;
  const Scalar _tmp205 =
      -_tmp197 * _tmp204 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp87), Scalar(2)) +
                     Scalar(0.77439999999999998) *
                         std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp84 - 1), Scalar(2))));
  const Scalar _tmp206 = Scalar(0.71007031138673404) * _tmp198;
  const Scalar _tmp207 = _tmp205 * _tmp206;
  const Scalar _tmp208 = Scalar(1.0) * _tmp203;
  const Scalar _tmp209 = _tmp194 + _tmp195 + _tmp196;
  const Scalar _tmp210 = std::pow(_tmp197, Scalar(-2));
  const Scalar _tmp211 = _tmp209 * _tmp210;
  const Scalar _tmp212 = Scalar(1.4083112389913199) * _tmp197;
  const Scalar _tmp213 =
      (_tmp198 * (_tmp135 * _tmp166 + _tmp154 * _tmp165 - _tmp199 - _tmp200 - _tmp201) -
       _tmp202 * _tmp211) /
      std::sqrt(Scalar(std::pow(_tmp202, Scalar(2)) * _tmp210 + 1));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = -Scalar(48.000000000000128) * _tmp0 -
               Scalar(1.4083112389913199) * fh1 *
                   (-Scalar(1.0) * _tmp38 * _tmp39 * fv1 * std::sinh(_tmp37) -
                    Scalar(34.083374946563197) * _tmp38 -
                    (-Scalar(0.71007031138673404) * _tmp34 * _tmp38 +
                     _tmp35 * (Scalar(1.4083112389913199) * _tmp31 * _tmp39 - _tmp33)) *
                        std::sinh(_tmp36)) +
               Scalar(1.4083112389913199) * std::cosh(_tmp36) -
               Scalar(1.4083112389913199) * std::cosh(_tmp37);
  _res(1, 0) =
      -_tmp159 *
          (-Scalar(34.083374946563197) * _tmp168 + Scalar(1.0) * _tmp169 * std::sinh(_tmp170) -
           (-Scalar(0.71007031138673404) * _tmp160 * _tmp168 +
            _tmp161 * (-_tmp158 * _tmp164 - _tmp159 * _tmp169)) *
               std::sinh(_tmp162)) -
      _tmp164 * (Scalar(34.083374946563197) * _tmp131 - std::cosh(_tmp162) + std::cosh(_tmp170));
  _res(2, 0) =
      -Scalar(1.4083112389913199) * _tmp174 *
          (Scalar(34.083374946563197) * _tmp177 - std::cosh(_tmp188) + std::cosh(_tmp189)) -
      _tmp193 *
          (-Scalar(34.083374946563197) * _tmp191 + Scalar(1.0) * _tmp192 * std::sinh(_tmp189) -
           (-Scalar(0.71007031138673404) * _tmp186 * _tmp191 +
            _tmp187 * (-_tmp174 * _tmp185 - _tmp192 * _tmp193)) *
               std::sinh(_tmp188));
  _res(3, 0) =
      -Scalar(1.4083112389913199) * _tmp209 *
          (Scalar(34.083374946563197) * _tmp198 - std::cosh(_tmp207) + std::cosh(_tmp208)) -
      _tmp212 *
          (-Scalar(34.083374946563197) * _tmp211 + Scalar(1.0) * _tmp213 * std::sinh(_tmp208) -
           (-Scalar(0.71007031138673404) * _tmp205 * _tmp211 +
            _tmp206 * (-_tmp204 * _tmp209 - _tmp212 * _tmp213)) *
               std::sinh(_tmp207));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
