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
 * Symbolic function: IK_residual_func_cost3_wrt_pa_Nl14
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
 *     res: Matrix43
 */
template <typename Scalar>
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost3WrtPaNl14(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 660

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (218)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp2 = 2 * _tmp1;
  const Scalar _tmp3 = _tmp0 * _tmp2;
  const Scalar _tmp4 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp5 = -_DeltaRot[0] * _Rot_init[0] - _DeltaRot[1] * _Rot_init[1] -
                       _DeltaRot[2] * _Rot_init[2] + _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp6 = 2 * _tmp5;
  const Scalar _tmp7 = _tmp4 * _tmp6;
  const Scalar _tmp8 = Scalar(0.20999999999999999) * _tmp3 - Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp9 = -_tmp8;
  const Scalar _tmp10 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp11 = 1 - 2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp12 = Scalar(0.20999999999999999) * _tmp10 + Scalar(0.20999999999999999) * _tmp11;
  const Scalar _tmp13 = _tmp2 * _tmp4;
  const Scalar _tmp14 = _tmp0 * _tmp6;
  const Scalar _tmp15 =
      -Scalar(0.010999999999999999) * _tmp13 - Scalar(0.010999999999999999) * _tmp14;
  const Scalar _tmp16 = -_tmp12 + _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp9;
  const Scalar _tmp18 = _tmp17 - p_a(0, 0) + position_vector(0, 0);
  const Scalar _tmp19 = 2 * _tmp0 * _tmp4;
  const Scalar _tmp20 = _tmp2 * _tmp5;
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp19 + Scalar(0.20999999999999999) * _tmp20;
  const Scalar _tmp22 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp23 =
      -Scalar(0.010999999999999999) * _tmp11 - Scalar(0.010999999999999999) * _tmp22;
  const Scalar _tmp24 = Scalar(0.20999999999999999) * _tmp13 - Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp25 = _tmp23 - _tmp24;
  const Scalar _tmp26 = _tmp21 + _tmp25;
  const Scalar _tmp27 = _tmp16 + _tmp8;
  const Scalar _tmp28 = _tmp27 - p_d(0, 0) + position_vector(0, 0);
  const Scalar _tmp29 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp7;
  const Scalar _tmp30 = -_tmp29;
  const Scalar _tmp31 = Scalar(0.20999999999999999) * _tmp10 +
                        Scalar(0.20999999999999999) * _tmp22 + Scalar(0.20999999999999999);
  const Scalar _tmp32 =
      -Scalar(0.010999999999999999) * _tmp19 + Scalar(0.010999999999999999) * _tmp20;
  const Scalar _tmp33 = _tmp31 + _tmp32;
  const Scalar _tmp34 = _tmp30 + _tmp33;
  const Scalar _tmp35 = _tmp34 - p_d(1, 0) + position_vector(1, 0);
  const Scalar _tmp36 = std::pow(Scalar(std::pow(_tmp28, Scalar(2)) + std::pow(_tmp35, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp37 = _tmp28 * _tmp36;
  const Scalar _tmp38 = -_tmp21;
  const Scalar _tmp39 = _tmp23 + _tmp24 + _tmp38;
  const Scalar _tmp40 = _tmp37 * _tmp39;
  const Scalar _tmp41 = _tmp35 * _tmp36;
  const Scalar _tmp42 = _tmp12 + _tmp15;
  const Scalar _tmp43 = _tmp42 + _tmp9;
  const Scalar _tmp44 = _tmp43 - p_b(0, 0) + position_vector(0, 0);
  const Scalar _tmp45 = Scalar(1.0) / (_tmp44);
  const Scalar _tmp46 = -_tmp31 + _tmp32;
  const Scalar _tmp47 = _tmp29 + _tmp46;
  const Scalar _tmp48 = _tmp47 - p_b(1, 0) + position_vector(1, 0);
  const Scalar _tmp49 = _tmp45 * _tmp48;
  const Scalar _tmp50 = std::pow(_tmp18, Scalar(2));
  const Scalar _tmp51 = _tmp30 + _tmp46;
  const Scalar _tmp52 = _tmp51 - p_a(1, 0) + position_vector(1, 0);
  const Scalar _tmp53 = std::pow(_tmp52, Scalar(2));
  const Scalar _tmp54 = _tmp50 + _tmp53;
  const Scalar _tmp55 = std::pow(_tmp54, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp56 = _tmp39 * _tmp55;
  const Scalar _tmp57 = _tmp18 * _tmp56;
  const Scalar _tmp58 = _tmp25 + _tmp38;
  const Scalar _tmp59 = _tmp52 * _tmp55;
  const Scalar _tmp60 = -_tmp49 * _tmp57 + _tmp58 * _tmp59;
  const Scalar _tmp61 = _tmp37 * _tmp49 - _tmp41;
  const Scalar _tmp62 = _tmp49 * _tmp55;
  const Scalar _tmp63 = _tmp18 * _tmp62 - _tmp59;
  const Scalar _tmp64 = Scalar(1.0) / (_tmp63);
  const Scalar _tmp65 = _tmp61 * _tmp64;
  const Scalar _tmp66 = Scalar(1.0) * _tmp47;
  const Scalar _tmp67 = Scalar(1.0) * _tmp43;
  const Scalar _tmp68 = (-_tmp17 + _tmp67) / (_tmp51 - _tmp66);
  const Scalar _tmp69 = _tmp55 * _tmp58;
  const Scalar _tmp70 = -_tmp18 * _tmp69 + _tmp57;
  const Scalar _tmp71 = -_tmp26 * _tmp37 + _tmp40 - _tmp65 * _tmp70 -
                        _tmp68 * (_tmp26 * _tmp41 - _tmp40 * _tmp49 - _tmp60 * _tmp65);
  const Scalar _tmp72 = Scalar(1.0) * _tmp64;
  const Scalar _tmp73 = _tmp68 * _tmp72;
  const Scalar _tmp74 = _tmp60 * _tmp73 - _tmp70 * _tmp72;
  const Scalar _tmp75 = Scalar(1.0) / (_tmp71);
  const Scalar _tmp76 = _tmp51 * _tmp55;
  const Scalar _tmp77 = _tmp17 * _tmp55;
  const Scalar _tmp78 =
      std::sqrt(Scalar(std::pow(_tmp44, Scalar(2)) + std::pow(_tmp48, Scalar(2))));
  const Scalar _tmp79 = Scalar(1.0) / (_tmp78);
  const Scalar _tmp80 = _tmp45 * _tmp78;
  const Scalar _tmp81 = _tmp80 * (_tmp43 * _tmp48 * _tmp79 - _tmp44 * _tmp47 * _tmp79);
  const Scalar _tmp82 = _tmp55 * _tmp81;
  const Scalar _tmp83 = _tmp18 * _tmp76 + _tmp18 * _tmp82 - _tmp52 * _tmp77;
  const Scalar _tmp84 = -_tmp27 * _tmp41 + _tmp34 * _tmp37 + _tmp37 * _tmp81 - _tmp65 * _tmp83;
  const Scalar _tmp85 = _tmp75 * _tmp84;
  const Scalar _tmp86 = -_tmp72 * _tmp83 - _tmp74 * _tmp85;
  const Scalar _tmp87 = Scalar(1.0) / (_tmp84);
  const Scalar _tmp88 = _tmp86 * _tmp87;
  const Scalar _tmp89 = _tmp71 * _tmp88 + _tmp74;
  const Scalar _tmp90 = _tmp61 * _tmp75;
  const Scalar _tmp91 = -_tmp89 * _tmp90 + Scalar(1.0);
  const Scalar _tmp92 = _tmp55 * _tmp64;
  const Scalar _tmp93 = _tmp91 * _tmp92;
  const Scalar _tmp94 = _tmp37 * _tmp75;
  const Scalar _tmp95 = _tmp42 + _tmp8;
  const Scalar _tmp96 = _tmp95 - p_c(0, 0) + position_vector(0, 0);
  const Scalar _tmp97 = _tmp29 + _tmp33;
  const Scalar _tmp98 = _tmp97 - p_c(1, 0) + position_vector(1, 0);
  const Scalar _tmp99 = std::pow(Scalar(std::pow(_tmp96, Scalar(2)) + std::pow(_tmp98, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp100 = _tmp98 * _tmp99;
  const Scalar _tmp101 = _tmp100 * fh1;
  const Scalar _tmp102 = _tmp101 * _tmp80;
  const Scalar _tmp103 = _tmp66 * _tmp68 + _tmp67;
  const Scalar _tmp104 = 0;
  const Scalar _tmp105 = _tmp104 * _tmp75;
  const Scalar _tmp106 = _tmp18 * _tmp55;
  const Scalar _tmp107 = _tmp105 * _tmp65;
  const Scalar _tmp108 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp109 = _tmp108 * _tmp80;
  const Scalar _tmp110 = _tmp49 * _tmp64;
  const Scalar _tmp111 = _tmp39 * _tmp49;
  const Scalar _tmp112 = _tmp110 * _tmp70 - _tmp39 - _tmp68 * (_tmp110 * _tmp60 + _tmp111);
  const Scalar _tmp113 = _tmp110 * _tmp83 - _tmp112 * _tmp85 - _tmp81;
  const Scalar _tmp114 = _tmp113 * _tmp87;
  const Scalar _tmp115 = _tmp112 + _tmp114 * _tmp71;
  const Scalar _tmp116 = -_tmp115 * _tmp90 - _tmp49;
  const Scalar _tmp117 = _tmp116 * _tmp92;
  const Scalar _tmp118 = _tmp96 * _tmp99;
  const Scalar _tmp119 = _tmp118 * fh1;
  const Scalar _tmp120 = _tmp119 * _tmp80;
  const Scalar _tmp121 = Scalar(1.0) * _tmp87;
  const Scalar _tmp122 = _tmp121 * _tmp65;
  const Scalar _tmp123 = fh1 * (_tmp100 * _tmp95 - _tmp118 * _tmp97);
  const Scalar _tmp124 = _tmp123 * _tmp80;
  const Scalar _tmp125 = std::exp(_tmp102 * (_tmp18 * _tmp93 + _tmp89 * _tmp94) +
                                  _tmp109 * (_tmp105 * _tmp37 - _tmp106 * _tmp107) +
                                  _tmp120 * (_tmp115 * _tmp94 + _tmp117 * _tmp18 + Scalar(1.0)) +
                                  _tmp124 * (-_tmp106 * _tmp122 + _tmp121 * _tmp37));
  const Scalar _tmp126 = std::pow(_tmp54, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp127 = _tmp126 * _tmp50;
  const Scalar _tmp128 = _tmp127 * _tmp64;
  const Scalar _tmp129 = _tmp127 * _tmp39;
  const Scalar _tmp130 = -_tmp127 * _tmp58 + _tmp129 - _tmp56 + _tmp69;
  const Scalar _tmp131 = _tmp126 * _tmp18 * _tmp52;
  const Scalar _tmp132 = _tmp127 * _tmp49 - _tmp131 - _tmp62;
  const Scalar _tmp133 = std::pow(_tmp63, Scalar(-2));
  const Scalar _tmp134 = _tmp133 * _tmp70;
  const Scalar _tmp135 = _tmp132 * _tmp134;
  const Scalar _tmp136 = _tmp133 * _tmp60;
  const Scalar _tmp137 = _tmp132 * _tmp61;
  const Scalar _tmp138 = _tmp131 * _tmp58;
  const Scalar _tmp139 = -_tmp129 * _tmp49 + _tmp138 + _tmp49 * _tmp56;
  const Scalar _tmp140 =
      -_tmp130 * _tmp65 + _tmp135 * _tmp61 - _tmp68 * (_tmp136 * _tmp137 - _tmp139 * _tmp65);
  const Scalar _tmp141 = std::pow(_tmp71, Scalar(-2));
  const Scalar _tmp142 = _tmp140 * _tmp141;
  const Scalar _tmp143 = _tmp142 * _tmp37;
  const Scalar _tmp144 = _tmp106 * _tmp91;
  const Scalar _tmp145 = _tmp132 * _tmp133;
  const Scalar _tmp146 = _tmp133 * _tmp83;
  const Scalar _tmp147 = _tmp127 * _tmp51 + _tmp127 * _tmp81 - _tmp131 * _tmp17 - _tmp76 - _tmp82;
  const Scalar _tmp148 = _tmp137 * _tmp146 - _tmp147 * _tmp65;
  const Scalar _tmp149 = std::pow(_tmp84, Scalar(-2));
  const Scalar _tmp150 = _tmp149 * _tmp71;
  const Scalar _tmp151 = _tmp150 * _tmp86;
  const Scalar _tmp152 = _tmp74 * _tmp75;
  const Scalar _tmp153 = _tmp132 * _tmp136;
  const Scalar _tmp154 =
      -_tmp130 * _tmp72 + Scalar(1.0) * _tmp135 + _tmp139 * _tmp73 - Scalar(1.0) * _tmp153 * _tmp68;
  const Scalar _tmp155 = _tmp142 * _tmp84;
  const Scalar _tmp156 = _tmp132 * _tmp146;
  const Scalar _tmp157 = _tmp71 * _tmp87;
  const Scalar _tmp158 = _tmp140 * _tmp88 - _tmp148 * _tmp151 + _tmp154 +
                         _tmp157 * (-_tmp147 * _tmp72 - _tmp148 * _tmp152 - _tmp154 * _tmp85 +
                                    _tmp155 * _tmp74 + Scalar(1.0) * _tmp156);
  const Scalar _tmp159 = _tmp61 * _tmp89;
  const Scalar _tmp160 = _tmp142 * _tmp159 - _tmp158 * _tmp90;
  const Scalar _tmp161 = _tmp18 * _tmp92;
  const Scalar _tmp162 = _tmp115 * _tmp61;
  const Scalar _tmp163 =
      _tmp110 * _tmp130 - _tmp135 * _tmp49 - _tmp68 * (_tmp110 * _tmp139 - _tmp153 * _tmp49);
  const Scalar _tmp164 = _tmp112 * _tmp75;
  const Scalar _tmp165 = _tmp113 * _tmp150;
  const Scalar _tmp166 = _tmp114 * _tmp140 - _tmp148 * _tmp165 +
                         _tmp157 * (_tmp110 * _tmp147 + _tmp112 * _tmp155 - _tmp148 * _tmp164 -
                                    _tmp156 * _tmp49 - _tmp163 * _tmp85) +
                         _tmp163;
  const Scalar _tmp167 = _tmp142 * _tmp162 - _tmp166 * _tmp90;
  const Scalar _tmp168 = _tmp116 * _tmp145;
  const Scalar _tmp169 = Scalar(1.0) * _tmp149;
  const Scalar _tmp170 = _tmp169 * _tmp37;
  const Scalar _tmp171 = _tmp133 * _tmp137;
  const Scalar _tmp172 = _tmp55 * _tmp65;
  const Scalar _tmp173 = _tmp106 * _tmp169 * _tmp65;
  const Scalar _tmp174 = _tmp104 * _tmp106;
  const Scalar _tmp175 = _tmp174 * _tmp65;
  const Scalar _tmp176 = _tmp174 * _tmp90;
  const Scalar _tmp177 = _tmp119 * _tmp64;
  const Scalar _tmp178 = _tmp101 * _tmp64;
  const Scalar _tmp179 = _tmp105 * _tmp108;
  const Scalar _tmp180 = _tmp121 * _tmp123;
  const Scalar _tmp181 =
      std::exp(-_tmp116 * _tmp177 - _tmp178 * _tmp91 + _tmp179 * _tmp65 + _tmp180 * _tmp65);
  const Scalar _tmp182 = _tmp104 * _tmp108;
  const Scalar _tmp183 = _tmp142 * _tmp182;
  const Scalar _tmp184 = _tmp101 * _tmp91;
  const Scalar _tmp185 = _tmp123 * _tmp169;
  const Scalar _tmp186 = _tmp148 * _tmp185;
  const Scalar _tmp187 = _tmp182 * _tmp90;
  const Scalar _tmp188 = _tmp101 * _tmp75;
  const Scalar _tmp189 = _tmp119 * _tmp75;
  const Scalar _tmp190 = std::exp(-_tmp115 * _tmp189 - _tmp179 - _tmp180 - _tmp188 * _tmp89);
  const Scalar _tmp191 = _tmp101 * _tmp89;
  const Scalar _tmp192 = _tmp115 * _tmp119;
  const Scalar _tmp193 = _tmp126 * _tmp53;
  const Scalar _tmp194 = _tmp131 * _tmp49 - _tmp193 + _tmp55;
  const Scalar _tmp195 = _tmp133 * _tmp194;
  const Scalar _tmp196 = _tmp116 * _tmp195;
  const Scalar _tmp197 = _tmp131 * _tmp64;
  const Scalar _tmp198 = _tmp194 * _tmp61;
  const Scalar _tmp199 = _tmp131 * _tmp39 - _tmp138;
  const Scalar _tmp200 = -_tmp111 * _tmp131 + _tmp193 * _tmp58 - _tmp69;
  const Scalar _tmp201 =
      _tmp134 * _tmp198 - _tmp199 * _tmp65 - _tmp68 * (_tmp136 * _tmp198 - _tmp200 * _tmp65);
  const Scalar _tmp202 = _tmp141 * _tmp201;
  const Scalar _tmp203 = _tmp202 * _tmp37;
  const Scalar _tmp204 = _tmp131 * _tmp51 + _tmp131 * _tmp81 - _tmp17 * _tmp193 + _tmp77;
  const Scalar _tmp205 = _tmp146 * _tmp198 - _tmp204 * _tmp65;
  const Scalar _tmp206 = _tmp194 * _tmp49;
  const Scalar _tmp207 =
      _tmp110 * _tmp199 - _tmp134 * _tmp206 - _tmp68 * (_tmp110 * _tmp200 - _tmp136 * _tmp206);
  const Scalar _tmp208 = _tmp202 * _tmp84;
  const Scalar _tmp209 = _tmp114 * _tmp201 +
                         _tmp157 * (_tmp110 * _tmp204 + _tmp112 * _tmp208 - _tmp146 * _tmp206 -
                                    _tmp164 * _tmp205 - _tmp207 * _tmp85) -
                         _tmp165 * _tmp205 + _tmp207;
  const Scalar _tmp210 = _tmp162 * _tmp202 - _tmp209 * _tmp90;
  const Scalar _tmp211 = Scalar(1.0) * _tmp194;
  const Scalar _tmp212 =
      _tmp134 * _tmp211 - _tmp136 * _tmp211 * _tmp68 - _tmp199 * _tmp72 + _tmp200 * _tmp73;
  const Scalar _tmp213 = -_tmp151 * _tmp205 +
                         _tmp157 * (_tmp146 * _tmp211 - _tmp152 * _tmp205 - _tmp204 * _tmp72 +
                                    _tmp208 * _tmp74 - _tmp212 * _tmp85) +
                         _tmp201 * _tmp88 + _tmp212;
  const Scalar _tmp214 = _tmp159 * _tmp202 - _tmp213 * _tmp90;
  const Scalar _tmp215 = _tmp133 * _tmp211 * _tmp61 * _tmp87;
  const Scalar _tmp216 = _tmp182 * _tmp202;
  const Scalar _tmp217 = _tmp185 * _tmp205;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) = -_tmp125 * (-_tmp102 * (_tmp128 * _tmp91 - _tmp143 * _tmp89 - _tmp144 * _tmp145 +
                                       _tmp158 * _tmp94 + _tmp160 * _tmp161 - _tmp93) -
                           _tmp109 * (-_tmp104 * _tmp143 + _tmp105 * _tmp172 - _tmp107 * _tmp127 +
                                      _tmp142 * _tmp175 + _tmp145 * _tmp176) -
                           _tmp120 * (-_tmp106 * _tmp168 - _tmp115 * _tmp143 + _tmp116 * _tmp128 -
                                      _tmp117 + _tmp161 * _tmp167 + _tmp166 * _tmp94) -
                           _tmp124 * (_tmp106 * _tmp121 * _tmp171 + _tmp121 * _tmp172 -
                                      _tmp122 * _tmp127 - _tmp148 * _tmp170 + _tmp148 * _tmp173));
  _res(2, 0) =
      -_tmp181 * (-_tmp119 * _tmp168 - _tmp145 * _tmp184 + _tmp145 * _tmp187 + _tmp160 * _tmp178 +
                  _tmp167 * _tmp177 + _tmp171 * _tmp180 + _tmp183 * _tmp65 + _tmp186 * _tmp65);
  _res(3, 0) = -_tmp190 * (-_tmp142 * _tmp191 - _tmp142 * _tmp192 + _tmp158 * _tmp188 +
                           _tmp166 * _tmp189 - _tmp183 - _tmp186);
  _res(0, 1) = 0;
  _res(1, 1) =
      -_tmp125 *
      (-_tmp102 * (-_tmp144 * _tmp195 + _tmp161 * _tmp214 + _tmp197 * _tmp91 - _tmp203 * _tmp89 +
                   _tmp213 * _tmp94) -
       _tmp109 * (-_tmp104 * _tmp203 - _tmp107 * _tmp131 + _tmp175 * _tmp202 + _tmp176 * _tmp195) -
       _tmp120 * (-_tmp106 * _tmp196 - _tmp115 * _tmp203 + _tmp116 * _tmp197 + _tmp161 * _tmp210 +
                  _tmp209 * _tmp94) -
       _tmp124 * (_tmp106 * _tmp215 - _tmp122 * _tmp131 - _tmp170 * _tmp205 + _tmp173 * _tmp205));
  _res(2, 1) =
      -_tmp181 * (-_tmp119 * _tmp196 + _tmp123 * _tmp215 + _tmp177 * _tmp210 + _tmp178 * _tmp214 -
                  _tmp184 * _tmp195 + _tmp187 * _tmp195 + _tmp216 * _tmp65 + _tmp217 * _tmp65);
  _res(3, 1) = -_tmp190 * (_tmp188 * _tmp213 + _tmp189 * _tmp209 - _tmp191 * _tmp202 -
                           _tmp192 * _tmp202 - _tmp216 - _tmp217);
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
