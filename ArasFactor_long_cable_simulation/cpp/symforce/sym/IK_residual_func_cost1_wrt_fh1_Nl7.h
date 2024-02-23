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
 * Symbolic function: IK_residual_func_cost1_wrt_fh1_Nl7
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1WrtFh1Nl7(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 653

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (212)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = _tmp0 * fv1;
  const Scalar _tmp2 = std::asinh(_tmp1);
  const Scalar _tmp3 = Scalar(1.0) * _tmp2;
  const Scalar _tmp4 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp5 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp6 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp7 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp8 = 1 - 2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp9 = Scalar(0.20999999999999999) * _tmp6 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp10 = -_tmp9;
  const Scalar _tmp11 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                        _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp12 = 2 * _tmp7;
  const Scalar _tmp13 = _tmp11 * _tmp12;
  const Scalar _tmp14 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                        2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp15 = _tmp14 * _tmp5;
  const Scalar _tmp16 = _tmp13 - _tmp15;
  const Scalar _tmp17 = -Scalar(0.010999999999999999) * _tmp16;
  const Scalar _tmp18 = 2 * _tmp11 * _tmp5;
  const Scalar _tmp19 = _tmp14 * _tmp7;
  const Scalar _tmp20 = Scalar(0.20999999999999999) * _tmp18 + Scalar(0.20999999999999999) * _tmp19;
  const Scalar _tmp21 = _tmp17 + _tmp20;
  const Scalar _tmp22 = _tmp10 + _tmp21;
  const Scalar _tmp23 = _tmp22 + position_vector(1, 0);
  const Scalar _tmp24 = -2 * std::pow(_tmp11, Scalar(2));
  const Scalar _tmp25 = Scalar(0.20999999999999999) * _tmp24 + Scalar(0.20999999999999999) * _tmp8;
  const Scalar _tmp26 = _tmp12 * _tmp5;
  const Scalar _tmp27 = _tmp11 * _tmp14;
  const Scalar _tmp28 = _tmp26 + _tmp27;
  const Scalar _tmp29 = -Scalar(0.010999999999999999) * _tmp28;
  const Scalar _tmp30 = Scalar(0.20999999999999999) * _tmp18 - Scalar(0.20999999999999999) * _tmp19;
  const Scalar _tmp31 = _tmp29 - _tmp30;
  const Scalar _tmp32 = _tmp25 + _tmp31;
  const Scalar _tmp33 = _tmp32 + position_vector(0, 0);
  const Scalar _tmp34 = Scalar(1.4083112389913199) * _tmp2;
  const Scalar _tmp35 =
      -_tmp34 * fh1 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp33), Scalar(2)) +
                     Scalar(0.77439999999999998) *
                         std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp23 - 1), Scalar(2))));
  const Scalar _tmp36 = Scalar(0.71007031138673404) * _tmp0;
  const Scalar _tmp37 = _tmp35 * _tmp36;
  const Scalar _tmp38 =
      std::pow(Scalar(_tmp4 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp39 = _tmp17 - _tmp20;
  const Scalar _tmp40 = _tmp10 + _tmp39;
  const Scalar _tmp41 = Scalar(1.0) * _tmp40;
  const Scalar _tmp42 = -_tmp41;
  const Scalar _tmp43 = _tmp39 + _tmp9;
  const Scalar _tmp44 = Scalar(1.0) / (_tmp42 + _tmp43);
  const Scalar _tmp45 = -_tmp25;
  const Scalar _tmp46 = _tmp29 + _tmp30;
  const Scalar _tmp47 = _tmp45 + _tmp46;
  const Scalar _tmp48 = _tmp31 + _tmp45;
  const Scalar _tmp49 = Scalar(1.0) * _tmp48;
  const Scalar _tmp50 = _tmp44 * (-_tmp47 + _tmp49);
  const Scalar _tmp51 = _tmp43 + position_vector(1, 0);
  const Scalar _tmp52 = _tmp51 + Scalar(-110.0);
  const Scalar _tmp53 = _tmp47 + position_vector(0, 0);
  const Scalar _tmp54 = _tmp53 + Scalar(125.0);
  const Scalar _tmp55 = std::pow(Scalar(std::pow(_tmp52, Scalar(2)) + std::pow(_tmp54, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp56 = _tmp52 * _tmp55;
  const Scalar _tmp57 = _tmp40 + position_vector(1, 0);
  const Scalar _tmp58 = _tmp57 + Scalar(110.0);
  const Scalar _tmp59 = _tmp48 + position_vector(0, 0);
  const Scalar _tmp60 = _tmp59 + Scalar(125.0);
  const Scalar _tmp61 = Scalar(1.0) / (_tmp60);
  const Scalar _tmp62 = _tmp58 * _tmp61;
  const Scalar _tmp63 = _tmp54 * _tmp55;
  const Scalar _tmp64 = Scalar(1.0) / (-_tmp56 + _tmp62 * _tmp63);
  const Scalar _tmp65 = Scalar(0.20999999999999999) * _tmp26 - Scalar(0.20999999999999999) * _tmp27;
  const Scalar _tmp66 = -_tmp65;
  const Scalar _tmp67 = -Scalar(0.010999999999999999) * _tmp24 -
                        Scalar(0.010999999999999999) * _tmp6 + Scalar(-0.010999999999999999);
  const Scalar _tmp68 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp69 = _tmp67 - _tmp68;
  const Scalar _tmp70 = _tmp66 + _tmp69;
  const Scalar _tmp71 = _tmp62 * _tmp70;
  const Scalar _tmp72 = _tmp67 + _tmp68;
  const Scalar _tmp73 = _tmp66 + _tmp72;
  const Scalar _tmp74 = _tmp64 * (_tmp56 * _tmp73 - _tmp63 * _tmp71);
  const Scalar _tmp75 = Scalar(1.0) * _tmp74;
  const Scalar _tmp76 = _tmp63 * _tmp70 - _tmp63 * _tmp73;
  const Scalar _tmp77 = Scalar(1.0) * _tmp64;
  const Scalar _tmp78 = _tmp50 * _tmp75 - _tmp76 * _tmp77;
  const Scalar _tmp79 = _tmp21 + _tmp9;
  const Scalar _tmp80 = _tmp79 + position_vector(1, 0);
  const Scalar _tmp81 = _tmp80 + Scalar(-110.0);
  const Scalar _tmp82 = _tmp25 + _tmp46;
  const Scalar _tmp83 = _tmp82 + position_vector(0, 0);
  const Scalar _tmp84 = _tmp83 + Scalar(-125.0);
  const Scalar _tmp85 = std::pow(Scalar(std::pow(_tmp81, Scalar(2)) + std::pow(_tmp84, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp86 = _tmp84 * _tmp85;
  const Scalar _tmp87 = _tmp65 + _tmp72;
  const Scalar _tmp88 = _tmp81 * _tmp85;
  const Scalar _tmp89 = _tmp62 * _tmp86 - _tmp88;
  const Scalar _tmp90 = _tmp64 * _tmp89;
  const Scalar _tmp91 = -_tmp71 * _tmp86 - _tmp74 * _tmp89 + _tmp87 * _tmp88;
  const Scalar _tmp92 = -_tmp50 * _tmp91 + _tmp70 * _tmp86 - _tmp76 * _tmp90 - _tmp86 * _tmp87;
  const Scalar _tmp93 = Scalar(1.0) / (_tmp92);
  const Scalar _tmp94 =
      std::sqrt(Scalar(std::pow(_tmp58, Scalar(2)) + std::pow(_tmp60, Scalar(2))));
  const Scalar _tmp95 = Scalar(1.0) / (_tmp94);
  const Scalar _tmp96 = _tmp61 * _tmp94;
  const Scalar _tmp97 = _tmp96 * (-_tmp40 * _tmp60 * _tmp95 + _tmp48 * _tmp58 * _tmp95);
  const Scalar _tmp98 = _tmp43 * _tmp63 - _tmp47 * _tmp56 + _tmp63 * _tmp97;
  const Scalar _tmp99 = _tmp79 * _tmp86 - _tmp82 * _tmp88 + _tmp86 * _tmp97 - _tmp90 * _tmp98;
  const Scalar _tmp100 = _tmp93 * _tmp99;
  const Scalar _tmp101 = Scalar(1.0) / (_tmp99);
  const Scalar _tmp102 = _tmp101 * _tmp92;
  const Scalar _tmp103 = _tmp102 * (-_tmp100 * _tmp78 - _tmp77 * _tmp98);
  const Scalar _tmp104 = _tmp93 * (_tmp103 + _tmp78);
  const Scalar _tmp105 = -_tmp104 * _tmp89 + Scalar(1.0);
  const Scalar _tmp106 = _tmp63 * _tmp64;
  const Scalar _tmp107 = _tmp23 + Scalar(110.0);
  const Scalar _tmp108 = _tmp33 + Scalar(-125.0);
  const Scalar _tmp109 =
      std::pow(Scalar(std::pow(_tmp107, Scalar(2)) + std::pow(_tmp108, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp110 = _tmp107 * _tmp109;
  const Scalar _tmp111 = _tmp110 * _tmp96 * (_tmp104 * _tmp86 + _tmp105 * _tmp106);
  const Scalar _tmp112 = _tmp108 * _tmp109;
  const Scalar _tmp113 = _tmp110 * _tmp32 - _tmp112 * _tmp22;
  const Scalar _tmp114 = Scalar(1.0) * _tmp101;
  const Scalar _tmp115 = _tmp63 * _tmp90;
  const Scalar _tmp116 = _tmp113 * _tmp96 * (-_tmp114 * _tmp115 + _tmp114 * _tmp86);
  const Scalar _tmp117 = _tmp62 * _tmp64;
  const Scalar _tmp118 = _tmp62 * _tmp74 + _tmp71;
  const Scalar _tmp119 = _tmp117 * _tmp76 - _tmp118 * _tmp50 - _tmp70;
  const Scalar _tmp120 = _tmp102 * (-_tmp100 * _tmp119 + _tmp117 * _tmp98 - _tmp97);
  const Scalar _tmp121 = _tmp93 * (_tmp119 + _tmp120);
  const Scalar _tmp122 = -_tmp121 * _tmp89 - _tmp62;
  const Scalar _tmp123 = _tmp112 * _tmp96 * (_tmp106 * _tmp122 + _tmp121 * _tmp86 + Scalar(1.0));
  const Scalar _tmp124 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp125 = _tmp41 * _tmp50 + _tmp49;
  const Scalar _tmp126 = 0;
  const Scalar _tmp127 = -_tmp111 * fh1 - _tmp116 * fh1 - _tmp123 * fh1 -
                         _tmp124 * _tmp96 * (-_tmp115 * _tmp126 + _tmp126 * _tmp86);
  const Scalar _tmp128 = std::pow(_tmp127, Scalar(-2));
  const Scalar _tmp129 = -_tmp111 - _tmp116 - _tmp123;
  const Scalar _tmp130 = _tmp128 * _tmp129;
  const Scalar _tmp131 = Scalar(1.0) / (_tmp127);
  const Scalar _tmp132 = _tmp42 + _tmp79;
  const Scalar _tmp133 = _tmp132 * _tmp50;
  const Scalar _tmp134 = Scalar(1.0) / (-_tmp133 + _tmp49 - _tmp82);
  const Scalar _tmp135 = _tmp125 * _tmp134;
  const Scalar _tmp136 = -_tmp126 * _tmp91 - _tmp132 * _tmp135 + _tmp42;
  const Scalar _tmp137 = Scalar(1.0) * _tmp44;
  const Scalar _tmp138 = Scalar(1.0) * _tmp134;
  const Scalar _tmp139 = _tmp132 * _tmp134;
  const Scalar _tmp140 = _tmp103 * _tmp139 - _tmp104 * _tmp91 - _tmp75;
  const Scalar _tmp141 = Scalar(1.0) * _tmp110 * (_tmp103 * _tmp138 - _tmp137 * _tmp140);
  const Scalar _tmp142 = _tmp65 + _tmp69;
  const Scalar _tmp143 = _tmp142 * fh1;
  const Scalar _tmp144 = -_tmp110 * _tmp143 - Scalar(40.024799999999999) * _tmp16 - _tmp22 * fv1;
  const Scalar _tmp145 = _tmp133 * _tmp138 + Scalar(1.0);
  const Scalar _tmp146 = _tmp138 * _tmp50;
  const Scalar _tmp147 = -Scalar(1.0) * _tmp137 * _tmp145 + Scalar(1.0) * _tmp146;
  const Scalar _tmp148 = _tmp112 * _tmp143 + Scalar(40.024799999999999) * _tmp28 + _tmp32 * fv1;
  const Scalar _tmp149 = _tmp132 * _tmp138;
  const Scalar _tmp150 = _tmp149 * _tmp44;
  const Scalar _tmp151 = -Scalar(1.0) * _tmp138 + Scalar(1.0) * _tmp150;
  const Scalar _tmp152 = _tmp102 * _tmp149 - _tmp114 * _tmp91;
  const Scalar _tmp153 = _tmp102 * _tmp138;
  const Scalar _tmp154 = Scalar(1.0) * _tmp113;
  const Scalar _tmp155 = _tmp154 * (-_tmp137 * _tmp152 + _tmp153);
  const Scalar _tmp156 = _tmp118 + _tmp120 * _tmp139 - _tmp121 * _tmp91;
  const Scalar _tmp157 = Scalar(1.0) * _tmp112 * (_tmp120 * _tmp138 - _tmp137 * _tmp156);
  const Scalar _tmp158 =
      Scalar(1.0) * _tmp124 * (-_tmp125 * _tmp138 - _tmp136 * _tmp137 + Scalar(1.0)) +
      _tmp141 * fh1 + _tmp144 * _tmp147 + _tmp148 * _tmp151 + _tmp155 * fh1 + _tmp157 * fh1;
  const Scalar _tmp159 = std::asinh(_tmp131 * _tmp158);
  const Scalar _tmp160 = Scalar(1.0) * _tmp159;
  const Scalar _tmp161 = _tmp110 * _tmp142;
  const Scalar _tmp162 = _tmp112 * _tmp142;
  const Scalar _tmp163 = (-_tmp130 * _tmp158 + _tmp131 * (_tmp141 - _tmp147 * _tmp161 +
                                                          _tmp151 * _tmp162 + _tmp155 + _tmp157)) /
                         std::sqrt(Scalar(_tmp128 * std::pow(_tmp158, Scalar(2)) + 1));
  const Scalar _tmp164 = Scalar(1.4083112389913199) * _tmp127;
  const Scalar _tmp165 =
      -_tmp159 * _tmp164 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(Scalar(0.77439999999999998) *
                         std::pow(Scalar(-Scalar(0.0090909090909090905) * _tmp57 - 1), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp59 - 1), Scalar(2))));
  const Scalar _tmp166 = Scalar(0.71007031138673404) * _tmp131;
  const Scalar _tmp167 = _tmp165 * _tmp166;
  const Scalar _tmp168 = Scalar(1.4083112389913199) * _tmp129;
  const Scalar _tmp169 = _tmp145 * _tmp44;
  const Scalar _tmp170 = _tmp113 * _tmp152 * _tmp44;
  const Scalar _tmp171 = _tmp110 * _tmp140 * _tmp44;
  const Scalar _tmp172 = _tmp112 * _tmp156 * _tmp44;
  const Scalar _tmp173 = _tmp138 * _tmp148;
  const Scalar _tmp174 = _tmp124 * _tmp136 * _tmp44 - _tmp132 * _tmp173 * _tmp44 +
                         _tmp144 * _tmp169 + _tmp170 * fh1 + _tmp171 * fh1 + _tmp172 * fh1;
  const Scalar _tmp175 = _tmp105 * _tmp110 * _tmp64;
  const Scalar _tmp176 = _tmp112 * _tmp122 * _tmp64;
  const Scalar _tmp177 = _tmp124 * _tmp126;
  const Scalar _tmp178 = _tmp101 * _tmp154;
  const Scalar _tmp179 = _tmp178 * fh1;
  const Scalar _tmp180 = _tmp175 * fh1 + _tmp176 * fh1 - _tmp177 * _tmp90 - _tmp179 * _tmp90;
  const Scalar _tmp181 = Scalar(1.0) / (_tmp180);
  const Scalar _tmp182 = std::asinh(_tmp174 * _tmp181);
  const Scalar _tmp183 = Scalar(1.0) * _tmp182;
  const Scalar _tmp184 = Scalar(1.4083112389913199) * _tmp180;
  const Scalar _tmp185 =
      -_tmp182 * _tmp184 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(Scalar(0.77439999999999998) *
                         std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp51), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.0080000000000000002) * _tmp53 - 1), Scalar(2))));
  const Scalar _tmp186 = Scalar(0.71007031138673404) * _tmp181;
  const Scalar _tmp187 = _tmp185 * _tmp186;
  const Scalar _tmp188 = _tmp175 + _tmp176 - _tmp178 * _tmp90;
  const Scalar _tmp189 = Scalar(1.4083112389913199) * _tmp188;
  const Scalar _tmp190 = std::pow(_tmp180, Scalar(-2));
  const Scalar _tmp191 = _tmp188 * _tmp190;
  const Scalar _tmp192 = (-_tmp174 * _tmp191 + _tmp181 * (-_tmp150 * _tmp162 - _tmp161 * _tmp169 +
                                                          _tmp170 + _tmp171 + _tmp172)) /
                         std::sqrt(Scalar(std::pow(_tmp174, Scalar(2)) * _tmp190 + 1));
  const Scalar _tmp193 = _tmp104 * _tmp110;
  const Scalar _tmp194 = _tmp112 * _tmp121;
  const Scalar _tmp195 = _tmp177 + _tmp179 + _tmp193 * fh1 + _tmp194 * fh1;
  const Scalar _tmp196 = Scalar(1.0) / (_tmp195);
  const Scalar _tmp197 = _tmp112 * _tmp120 * _tmp134;
  const Scalar _tmp198 = _tmp103 * _tmp110 * _tmp134;
  const Scalar _tmp199 = _tmp113 * _tmp153;
  const Scalar _tmp200 = _tmp124 * _tmp135 - _tmp144 * _tmp146 + _tmp173 - _tmp197 * fh1 -
                         _tmp198 * fh1 - _tmp199 * fh1;
  const Scalar _tmp201 = std::asinh(_tmp196 * _tmp200);
  const Scalar _tmp202 = Scalar(1.4083112389913199) * _tmp195;
  const Scalar _tmp203 =
      -_tmp201 * _tmp202 -
      Scalar(125.0) *
          std::sqrt(
              Scalar(Scalar(0.77439999999999998) *
                         std::pow(Scalar(1 - Scalar(0.0090909090909090905) * _tmp80), Scalar(2)) +
                     std::pow(Scalar(1 - Scalar(0.0080000000000000002) * _tmp83), Scalar(2))));
  const Scalar _tmp204 = Scalar(0.71007031138673404) * _tmp196;
  const Scalar _tmp205 = _tmp203 * _tmp204;
  const Scalar _tmp206 = Scalar(1.0) * _tmp201;
  const Scalar _tmp207 = _tmp178 + _tmp193 + _tmp194;
  const Scalar _tmp208 = Scalar(1.4083112389913199) * _tmp207;
  const Scalar _tmp209 = std::pow(_tmp195, Scalar(-2));
  const Scalar _tmp210 = _tmp207 * _tmp209;
  const Scalar _tmp211 =
      (_tmp196 * (_tmp138 * _tmp162 + _tmp146 * _tmp161 - _tmp197 - _tmp198 - _tmp199) -
       _tmp200 * _tmp210) /
      std::sqrt(Scalar(std::pow(_tmp200, Scalar(2)) * _tmp209 + 1));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = -Scalar(48.000000000000128) * _tmp0 -
               Scalar(1.4083112389913199) * fh1 *
                   (-Scalar(1.0) * _tmp38 * _tmp4 * fv1 * std::sinh(_tmp3) -
                    Scalar(34.083374946563197) * _tmp4 -
                    (-Scalar(0.71007031138673404) * _tmp35 * _tmp4 +
                     _tmp36 * (Scalar(1.4083112389913199) * _tmp1 * _tmp38 - _tmp34)) *
                        std::sinh(_tmp37)) -
               Scalar(1.4083112389913199) * std::cosh(_tmp3) +
               Scalar(1.4083112389913199) * std::cosh(_tmp37);
  _res(1, 0) =
      -_tmp164 *
          (-Scalar(34.083374946563197) * _tmp130 + Scalar(1.0) * _tmp163 * std::sinh(_tmp160) -
           (-Scalar(0.71007031138673404) * _tmp130 * _tmp165 +
            _tmp166 * (-_tmp159 * _tmp168 - _tmp163 * _tmp164)) *
               std::sinh(_tmp167)) -
      _tmp168 * (Scalar(34.083374946563197) * _tmp131 + std::cosh(_tmp160) - std::cosh(_tmp167));
  _res(2, 0) =
      -_tmp184 *
          (-Scalar(34.083374946563197) * _tmp191 + Scalar(1.0) * _tmp192 * std::sinh(_tmp183) -
           (-Scalar(0.71007031138673404) * _tmp185 * _tmp191 +
            _tmp186 * (-_tmp182 * _tmp189 - _tmp184 * _tmp192)) *
               std::sinh(_tmp187)) -
      _tmp189 * (Scalar(34.083374946563197) * _tmp181 + std::cosh(_tmp183) - std::cosh(_tmp187));
  _res(3, 0) =
      -_tmp202 *
          (-Scalar(34.083374946563197) * _tmp210 + Scalar(1.0) * _tmp211 * std::sinh(_tmp206) -
           (-Scalar(0.71007031138673404) * _tmp203 * _tmp210 +
            _tmp204 * (-_tmp201 * _tmp208 - _tmp202 * _tmp211)) *
               std::sinh(_tmp205)) -
      _tmp208 * (Scalar(34.083374946563197) * _tmp196 - std::cosh(_tmp205) + std::cosh(_tmp206));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
