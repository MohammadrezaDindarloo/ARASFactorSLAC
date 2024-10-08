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
 * Symbolic function: IK_residual_func_cost3_wrt_pd_Nl17
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost3WrtPdNl17(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 723

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (245)
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
  const Scalar _tmp7 = Scalar(0.20999999999999999) * _tmp3 - Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp8 = 2 * _tmp0 * _tmp4;
  const Scalar _tmp9 = _tmp1 * _tmp5;
  const Scalar _tmp10 =
      -Scalar(0.010999999999999999) * _tmp8 - Scalar(0.010999999999999999) * _tmp9;
  const Scalar _tmp11 = -2 * std::pow(_tmp4, Scalar(2));
  const Scalar _tmp12 = 1 - 2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp13 = Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp14 = _tmp10 - _tmp13;
  const Scalar _tmp15 = _tmp14 + _tmp7;
  const Scalar _tmp16 = _tmp15 - p_d(0, 0) + position_vector(0, 0);
  const Scalar _tmp17 = Scalar(1.0) / (_tmp16);
  const Scalar _tmp18 = std::pow(_tmp16, Scalar(2));
  const Scalar _tmp19 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp20 = Scalar(0.20999999999999999) * _tmp11 +
                        Scalar(0.20999999999999999) * _tmp19 + Scalar(0.20999999999999999);
  const Scalar _tmp21 = _tmp2 * _tmp4;
  const Scalar _tmp22 = _tmp0 * _tmp5;
  const Scalar _tmp23 =
      -Scalar(0.010999999999999999) * _tmp21 + Scalar(0.010999999999999999) * _tmp22;
  const Scalar _tmp24 = Scalar(0.20999999999999999) * _tmp3 + Scalar(0.20999999999999999) * _tmp6;
  const Scalar _tmp25 = _tmp23 - _tmp24;
  const Scalar _tmp26 = _tmp20 + _tmp25;
  const Scalar _tmp27 = _tmp26 - p_d(1, 0) + position_vector(1, 0);
  const Scalar _tmp28 = std::pow(_tmp27, Scalar(2));
  const Scalar _tmp29 = _tmp18 + _tmp28;
  const Scalar _tmp30 = std::sqrt(_tmp29);
  const Scalar _tmp31 = _tmp17 * _tmp30;
  const Scalar _tmp32 = Scalar(0.20999999999999999) * _tmp8 - Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp33 =
      -Scalar(0.010999999999999999) * _tmp12 - Scalar(0.010999999999999999) * _tmp19;
  const Scalar _tmp34 = Scalar(0.20999999999999999) * _tmp21 + Scalar(0.20999999999999999) * _tmp22;
  const Scalar _tmp35 = _tmp33 - _tmp34;
  const Scalar _tmp36 = _tmp32 + _tmp35;
  const Scalar _tmp37 = -_tmp7;
  const Scalar _tmp38 = _tmp10 + _tmp13;
  const Scalar _tmp39 = _tmp37 + _tmp38;
  const Scalar _tmp40 = _tmp39 - p_b(0, 0) + position_vector(0, 0);
  const Scalar _tmp41 = std::pow(_tmp40, Scalar(2));
  const Scalar _tmp42 = -_tmp20;
  const Scalar _tmp43 = _tmp23 + _tmp24;
  const Scalar _tmp44 = _tmp42 + _tmp43;
  const Scalar _tmp45 = _tmp44 - p_b(1, 0) + position_vector(1, 0);
  const Scalar _tmp46 = _tmp41 + std::pow(_tmp45, Scalar(2));
  const Scalar _tmp47 = std::pow(_tmp46, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp48 = _tmp45 * _tmp47;
  const Scalar _tmp49 = _tmp17 * _tmp27;
  const Scalar _tmp50 = -_tmp32;
  const Scalar _tmp51 = _tmp33 + _tmp34 + _tmp50;
  const Scalar _tmp52 = _tmp40 * _tmp47;
  const Scalar _tmp53 = _tmp51 * _tmp52;
  const Scalar _tmp54 = _tmp36 * _tmp48 - _tmp49 * _tmp53;
  const Scalar _tmp55 = Scalar(1.0) * _tmp15;
  const Scalar _tmp56 = Scalar(1.0) * _tmp26;
  const Scalar _tmp57 = (-_tmp39 + _tmp55) / (_tmp44 - _tmp56);
  const Scalar _tmp58 = -_tmp48 + _tmp49 * _tmp52;
  const Scalar _tmp59 = Scalar(1.0) / (_tmp58);
  const Scalar _tmp60 = Scalar(1.0) * _tmp59;
  const Scalar _tmp61 = _tmp57 * _tmp60;
  const Scalar _tmp62 = -_tmp36 * _tmp52 + _tmp53;
  const Scalar _tmp63 = _tmp54 * _tmp61 - _tmp60 * _tmp62;
  const Scalar _tmp64 = _tmp14 + _tmp37;
  const Scalar _tmp65 = _tmp64 - p_a(0, 0) + position_vector(0, 0);
  const Scalar _tmp66 = _tmp25 + _tmp42;
  const Scalar _tmp67 = _tmp66 - p_a(1, 0) + position_vector(1, 0);
  const Scalar _tmp68 = std::pow(Scalar(std::pow(_tmp65, Scalar(2)) + std::pow(_tmp67, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp69 = _tmp65 * _tmp68;
  const Scalar _tmp70 = Scalar(1.0) / (_tmp30);
  const Scalar _tmp71 = _tmp15 * _tmp70;
  const Scalar _tmp72 = _tmp26 * _tmp70;
  const Scalar _tmp73 = -_tmp16 * _tmp72 + _tmp27 * _tmp71;
  const Scalar _tmp74 = _tmp30 * _tmp73;
  const Scalar _tmp75 = _tmp17 * _tmp74;
  const Scalar _tmp76 = -_tmp39 * _tmp48 + _tmp44 * _tmp52 + _tmp52 * _tmp75;
  const Scalar _tmp77 = _tmp67 * _tmp68;
  const Scalar _tmp78 = _tmp49 * _tmp69 - _tmp77;
  const Scalar _tmp79 = _tmp59 * _tmp78;
  const Scalar _tmp80 = -_tmp64 * _tmp77 + _tmp66 * _tmp69 + _tmp69 * _tmp75 - _tmp76 * _tmp79;
  const Scalar _tmp81 = _tmp35 + _tmp50;
  const Scalar _tmp82 = _tmp49 * _tmp51;
  const Scalar _tmp83 = _tmp51 * _tmp69 -
                        _tmp57 * (-_tmp54 * _tmp79 - _tmp69 * _tmp82 + _tmp77 * _tmp81) -
                        _tmp62 * _tmp79 - _tmp69 * _tmp81;
  const Scalar _tmp84 = Scalar(1.0) / (_tmp83);
  const Scalar _tmp85 = _tmp80 * _tmp84;
  const Scalar _tmp86 = -_tmp60 * _tmp76 - _tmp63 * _tmp85;
  const Scalar _tmp87 = Scalar(1.0) / (_tmp80);
  const Scalar _tmp88 = _tmp83 * _tmp87;
  const Scalar _tmp89 = _tmp63 + _tmp86 * _tmp88;
  const Scalar _tmp90 = _tmp78 * _tmp84;
  const Scalar _tmp91 = -_tmp89 * _tmp90 + Scalar(1.0);
  const Scalar _tmp92 = _tmp52 * _tmp59;
  const Scalar _tmp93 = _tmp69 * _tmp84;
  const Scalar _tmp94 = _tmp89 * _tmp93;
  const Scalar _tmp95 = _tmp38 + _tmp7;
  const Scalar _tmp96 = _tmp95 - p_c(0, 0) + position_vector(0, 0);
  const Scalar _tmp97 = _tmp20 + _tmp43;
  const Scalar _tmp98 = _tmp97 - p_c(1, 0) + position_vector(1, 0);
  const Scalar _tmp99 = std::pow(Scalar(std::pow(_tmp96, Scalar(2)) + std::pow(_tmp98, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp100 = _tmp98 * _tmp99;
  const Scalar _tmp101 = _tmp100 * (_tmp91 * _tmp92 + _tmp94);
  const Scalar _tmp102 = _tmp101 * fh1;
  const Scalar _tmp103 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp104 = _tmp55 + _tmp56 * _tmp57;
  const Scalar _tmp105 = 0;
  const Scalar _tmp106 = _tmp105 * _tmp84;
  const Scalar _tmp107 = _tmp52 * _tmp79;
  const Scalar _tmp108 = _tmp106 * _tmp69;
  const Scalar _tmp109 = _tmp103 * (-_tmp106 * _tmp107 + _tmp108);
  const Scalar _tmp110 = _tmp49 * _tmp59;
  const Scalar _tmp111 = _tmp110 * _tmp62 - _tmp51 - _tmp57 * (_tmp110 * _tmp54 + _tmp82);
  const Scalar _tmp112 = _tmp110 * _tmp76 - _tmp111 * _tmp85 - _tmp75;
  const Scalar _tmp113 = _tmp111 + _tmp112 * _tmp88;
  const Scalar _tmp114 = -_tmp113 * _tmp90 - _tmp49;
  const Scalar _tmp115 = _tmp113 * _tmp93;
  const Scalar _tmp116 = _tmp96 * _tmp99;
  const Scalar _tmp117 = _tmp116 * (_tmp114 * _tmp92 + _tmp115 + Scalar(1.0));
  const Scalar _tmp118 = _tmp117 * fh1;
  const Scalar _tmp119 = Scalar(1.0) * _tmp87;
  const Scalar _tmp120 = fh1 * (_tmp100 * _tmp95 - _tmp116 * _tmp97);
  const Scalar _tmp121 = _tmp120 * (-_tmp107 * _tmp119 + _tmp119 * _tmp69);
  const Scalar _tmp122 =
      std::exp(_tmp102 * _tmp31 + _tmp109 * _tmp31 + _tmp118 * _tmp31 + _tmp121 * _tmp31);
  const Scalar _tmp123 = _tmp70 * fh1;
  const Scalar _tmp124 = _tmp101 * _tmp123;
  const Scalar _tmp125 = _tmp109 * _tmp70;
  const Scalar _tmp126 = Scalar(1.0) / (_tmp18);
  const Scalar _tmp127 = _tmp126 * _tmp30;
  const Scalar _tmp128 = _tmp126 * _tmp27;
  const Scalar _tmp129 = _tmp60 * _tmp69 * _tmp87;
  const Scalar _tmp130 = _tmp128 * _tmp129;
  const Scalar _tmp131 = _tmp41 / _tmp46;
  const Scalar _tmp132 = std::pow(_tmp58, Scalar(-2));
  const Scalar _tmp133 = _tmp128 * _tmp132;
  const Scalar _tmp134 = _tmp131 * _tmp133;
  const Scalar _tmp135 = _tmp134 * _tmp78;
  const Scalar _tmp136 = _tmp70 * _tmp73;
  const Scalar _tmp137 = _tmp128 * _tmp59;
  const Scalar _tmp138 = _tmp137 * _tmp76;
  const Scalar _tmp139 = _tmp126 * _tmp74;
  const Scalar _tmp140 = std::pow(_tmp29, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp141 = _tmp140 * _tmp26;
  const Scalar _tmp142 = _tmp140 * _tmp15;
  const Scalar _tmp143 = _tmp16 * _tmp27;
  const Scalar _tmp144 = _tmp31 * (-_tmp141 * _tmp18 + _tmp142 * _tmp143 + _tmp72);
  const Scalar _tmp145 = -_tmp136 * _tmp52 + _tmp139 * _tmp52 + _tmp144 * _tmp52;
  const Scalar _tmp146 = _tmp52 * _tmp76;
  const Scalar _tmp147 = _tmp133 * _tmp146;
  const Scalar _tmp148 = -_tmp136 * _tmp69 - _tmp138 * _tmp69 + _tmp139 * _tmp69 +
                         _tmp144 * _tmp69 - _tmp145 * _tmp79 + _tmp147 * _tmp78;
  const Scalar _tmp149 = std::pow(_tmp80, Scalar(-2));
  const Scalar _tmp150 = Scalar(1.0) * _tmp149;
  const Scalar _tmp151 = _tmp150 * _tmp69;
  const Scalar _tmp152 = _tmp107 * _tmp150;
  const Scalar _tmp153 = _tmp120 * _tmp31;
  const Scalar _tmp154 = _tmp128 * _tmp53;
  const Scalar _tmp155 = _tmp52 * _tmp54;
  const Scalar _tmp156 = _tmp133 * _tmp155;
  const Scalar _tmp157 = _tmp128 * _tmp51;
  const Scalar _tmp158 = _tmp137 * _tmp54;
  const Scalar _tmp159 = _tmp137 * _tmp62;
  const Scalar _tmp160 = _tmp52 * _tmp62;
  const Scalar _tmp161 = _tmp132 * _tmp160;
  const Scalar _tmp162 = _tmp128 * _tmp161;
  const Scalar _tmp163 =
      -_tmp159 * _tmp69 + _tmp162 * _tmp78 -
      _tmp57 * (_tmp154 * _tmp79 + _tmp156 * _tmp78 - _tmp157 * _tmp69 - _tmp158 * _tmp69);
  const Scalar _tmp164 = std::pow(_tmp83, Scalar(-2));
  const Scalar _tmp165 = _tmp164 * _tmp78;
  const Scalar _tmp166 = _tmp163 * _tmp165;
  const Scalar _tmp167 = _tmp149 * _tmp83;
  const Scalar _tmp168 = _tmp148 * _tmp167;
  const Scalar _tmp169 = _tmp164 * _tmp80;
  const Scalar _tmp170 = _tmp111 * _tmp169;
  const Scalar _tmp171 = _tmp111 * _tmp84;
  const Scalar _tmp172 = _tmp28 / [&]() {
    const Scalar base = _tmp16;
    return base * base * base;
  }();
  const Scalar _tmp173 = _tmp132 * _tmp172;
  const Scalar _tmp174 =
      _tmp159 - _tmp161 * _tmp172 -
      _tmp57 * (-_tmp155 * _tmp173 + _tmp157 + _tmp158 - _tmp172 * _tmp53 * _tmp59);
  const Scalar _tmp175 = _tmp163 * _tmp87;
  const Scalar _tmp176 =
      -_tmp112 * _tmp168 + _tmp112 * _tmp175 + _tmp174 +
      _tmp88 * (_tmp110 * _tmp145 + _tmp136 + _tmp138 - _tmp139 - _tmp144 - _tmp146 * _tmp173 -
                _tmp148 * _tmp171 + _tmp163 * _tmp170 - _tmp174 * _tmp85);
  const Scalar _tmp177 = _tmp113 * _tmp166 - _tmp115 * _tmp128 - _tmp128 - _tmp176 * _tmp90;
  const Scalar _tmp178 = _tmp163 * _tmp164;
  const Scalar _tmp179 = _tmp113 * _tmp69;
  const Scalar _tmp180 = _tmp116 * fh1;
  const Scalar _tmp181 = _tmp180 * _tmp31;
  const Scalar _tmp182 = _tmp121 * _tmp70;
  const Scalar _tmp183 = _tmp117 * _tmp123;
  const Scalar _tmp184 = _tmp178 * _tmp69;
  const Scalar _tmp185 = _tmp105 * _tmp107;
  const Scalar _tmp186 = _tmp103 * _tmp31;
  const Scalar _tmp187 = _tmp169 * _tmp63;
  const Scalar _tmp188 = _tmp63 * _tmp84;
  const Scalar _tmp189 = -_tmp154 * _tmp61 - Scalar(1.0) * _tmp156 * _tmp57 + Scalar(1.0) * _tmp162;
  const Scalar _tmp190 = -_tmp168 * _tmp86 + _tmp175 * _tmp86 + _tmp189 +
                         _tmp88 * (-_tmp145 * _tmp60 + Scalar(1.0) * _tmp147 - _tmp148 * _tmp188 +
                                   _tmp163 * _tmp187 - _tmp189 * _tmp85);
  const Scalar _tmp191 = _tmp131 * _tmp91;
  const Scalar _tmp192 = -_tmp128 * _tmp94 + _tmp166 * _tmp89 - _tmp190 * _tmp90;
  const Scalar _tmp193 = _tmp100 * fh1;
  const Scalar _tmp194 = _tmp193 * _tmp31;
  const Scalar _tmp195 = _tmp59 * fh1;
  const Scalar _tmp196 = _tmp116 * _tmp195;
  const Scalar _tmp197 = _tmp103 * _tmp106;
  const Scalar _tmp198 = _tmp119 * _tmp120;
  const Scalar _tmp199 = _tmp100 * _tmp195;
  const Scalar _tmp200 =
      std::exp(-_tmp114 * _tmp196 + _tmp197 * _tmp79 + _tmp198 * _tmp79 - _tmp199 * _tmp91);
  const Scalar _tmp201 = _tmp197 * _tmp69;
  const Scalar _tmp202 = _tmp193 * _tmp52 * _tmp91;
  const Scalar _tmp203 = _tmp120 * _tmp150;
  const Scalar _tmp204 = _tmp148 * _tmp203;
  const Scalar _tmp205 = _tmp180 * _tmp52;
  const Scalar _tmp206 = _tmp133 * _tmp52 * _tmp78;
  const Scalar _tmp207 = _tmp103 * _tmp105;
  const Scalar _tmp208 = _tmp178 * _tmp207;
  const Scalar _tmp209 = _tmp84 * fh1;
  const Scalar _tmp210 = _tmp100 * _tmp209;
  const Scalar _tmp211 = _tmp116 * _tmp209;
  const Scalar _tmp212 = std::exp(-_tmp113 * _tmp211 - _tmp197 - _tmp198 - _tmp210 * _tmp89);
  const Scalar _tmp213 = _tmp113 * _tmp180;
  const Scalar _tmp214 = _tmp193 * _tmp89;
  const Scalar _tmp215 = _tmp136 * _tmp49;
  const Scalar _tmp216 = _tmp132 * _tmp17;
  const Scalar _tmp217 = _tmp216 * _tmp78;
  const Scalar _tmp218 = _tmp31 * (-_tmp141 * _tmp143 + _tmp142 * _tmp28 - _tmp71);
  const Scalar _tmp219 = -_tmp215 * _tmp52 + _tmp218 * _tmp52;
  const Scalar _tmp220 = _tmp17 * _tmp59;
  const Scalar _tmp221 = _tmp220 * _tmp76;
  const Scalar _tmp222 = -_tmp146 * _tmp217 - _tmp215 * _tmp69 + _tmp218 * _tmp69 -
                         _tmp219 * _tmp79 + _tmp221 * _tmp69;
  const Scalar _tmp223 = _tmp220 * _tmp62;
  const Scalar _tmp224 = _tmp17 * _tmp51;
  const Scalar _tmp225 = _tmp220 * _tmp54;
  const Scalar _tmp226 =
      _tmp162 - _tmp223 - _tmp57 * (_tmp137 * _tmp53 + _tmp156 - _tmp224 - _tmp225);
  const Scalar _tmp227 =
      -_tmp160 * _tmp217 + _tmp223 * _tmp69 -
      _tmp57 * (-_tmp107 * _tmp224 - _tmp155 * _tmp217 + _tmp224 * _tmp69 + _tmp225 * _tmp69);
  const Scalar _tmp228 = _tmp227 * _tmp87;
  const Scalar _tmp229 = _tmp167 * _tmp222;
  const Scalar _tmp230 =
      _tmp112 * _tmp228 - _tmp112 * _tmp229 + _tmp226 +
      _tmp88 * (_tmp110 * _tmp219 + _tmp147 + _tmp170 * _tmp227 - _tmp171 * _tmp222 + _tmp215 -
                _tmp218 - _tmp221 - _tmp226 * _tmp85);
  const Scalar _tmp231 = _tmp165 * _tmp227;
  const Scalar _tmp232 = _tmp113 * _tmp231 + _tmp115 * _tmp17 + _tmp17 - _tmp230 * _tmp90;
  const Scalar _tmp233 = _tmp114 * _tmp216;
  const Scalar _tmp234 = _tmp164 * _tmp227;
  const Scalar _tmp235 = _tmp129 * _tmp17;
  const Scalar _tmp236 = _tmp131 * _tmp217;
  const Scalar _tmp237 = _tmp234 * _tmp69;
  const Scalar _tmp238 = Scalar(1.0) * _tmp216;
  const Scalar _tmp239 = _tmp155 * _tmp238 * _tmp57 - _tmp160 * _tmp238 + _tmp224 * _tmp52 * _tmp61;
  const Scalar _tmp240 = _tmp228 * _tmp86 - _tmp229 * _tmp86 + _tmp239 +
                         _tmp88 * (-_tmp146 * _tmp238 + _tmp187 * _tmp227 - _tmp188 * _tmp222 -
                                   _tmp219 * _tmp60 - _tmp239 * _tmp85);
  const Scalar _tmp241 = _tmp17 * _tmp94 + _tmp231 * _tmp89 - _tmp240 * _tmp90;
  const Scalar _tmp242 = _tmp217 * _tmp52;
  const Scalar _tmp243 = _tmp207 * _tmp234;
  const Scalar _tmp244 = _tmp203 * _tmp222;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) =
      -_tmp122 *
      (-_tmp102 * _tmp127 - _tmp109 * _tmp127 - _tmp118 * _tmp127 - _tmp121 * _tmp127 + _tmp124 +
       _tmp125 -
       _tmp153 * (_tmp119 * _tmp135 - _tmp130 * _tmp52 - _tmp148 * _tmp151 + _tmp148 * _tmp152) -
       _tmp181 * (-_tmp114 * _tmp134 + _tmp176 * _tmp93 + _tmp177 * _tmp92 - _tmp178 * _tmp179) +
       _tmp182 + _tmp183 -
       _tmp186 * (-_tmp105 * _tmp184 + _tmp106 * _tmp135 - _tmp108 * _tmp128 * _tmp92 +
                  _tmp178 * _tmp185) -
       _tmp194 * (-_tmp133 * _tmp191 - _tmp184 * _tmp89 + _tmp190 * _tmp93 + _tmp192 * _tmp92));
  _res(2, 0) =
      -_tmp200 * (-_tmp114 * _tmp133 * _tmp205 - _tmp120 * _tmp130 - _tmp133 * _tmp202 -
                  _tmp137 * _tmp201 + _tmp177 * _tmp196 + _tmp192 * _tmp199 + _tmp197 * _tmp206 +
                  _tmp198 * _tmp206 + _tmp204 * _tmp79 + _tmp208 * _tmp79);
  _res(3, 0) = -_tmp212 * (_tmp176 * _tmp211 - _tmp178 * _tmp213 - _tmp178 * _tmp214 +
                           _tmp190 * _tmp210 - _tmp204 - _tmp208);
  _res(0, 1) = 0;
  _res(1, 1) =
      -_tmp122 *
      (_tmp124 * _tmp49 + _tmp125 * _tmp49 -
       _tmp153 * (-_tmp119 * _tmp236 - _tmp151 * _tmp222 + _tmp152 * _tmp222 + _tmp235 * _tmp52) -
       _tmp181 * (_tmp131 * _tmp233 - _tmp179 * _tmp234 + _tmp230 * _tmp93 + _tmp232 * _tmp92) +
       _tmp182 * _tmp49 + _tmp183 * _tmp49 -
       _tmp186 * (-_tmp105 * _tmp237 - _tmp106 * _tmp236 + _tmp108 * _tmp220 * _tmp52 +
                  _tmp185 * _tmp234) -
       _tmp194 * (_tmp191 * _tmp216 - _tmp237 * _tmp89 + _tmp240 * _tmp93 + _tmp241 * _tmp92));
  _res(2, 1) =
      -_tmp200 * (_tmp120 * _tmp235 + _tmp196 * _tmp232 - _tmp197 * _tmp242 - _tmp198 * _tmp242 +
                  _tmp199 * _tmp241 + _tmp201 * _tmp220 + _tmp202 * _tmp216 + _tmp205 * _tmp233 +
                  _tmp243 * _tmp79 + _tmp244 * _tmp79);
  _res(3, 1) = -_tmp212 * (_tmp210 * _tmp240 + _tmp211 * _tmp230 - _tmp213 * _tmp234 -
                           _tmp214 * _tmp234 - _tmp243 - _tmp244);
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
