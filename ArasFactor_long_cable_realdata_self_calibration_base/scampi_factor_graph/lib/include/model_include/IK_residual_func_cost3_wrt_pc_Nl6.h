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
 * Symbolic function: IK_residual_func_cost3_wrt_pc_Nl6
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
Eigen::Matrix<Scalar, 4, 3> IkResidualFuncCost3WrtPcNl6(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const Eigen::Matrix<Scalar, 4, 1>& encoder,
    const Eigen::Matrix<Scalar, 3, 1>& p_a, const Eigen::Matrix<Scalar, 3, 1>& p_b,
    const Eigen::Matrix<Scalar, 3, 1>& p_c, const Eigen::Matrix<Scalar, 3, 1>& p_d,
    const sym::Rot3<Scalar>& Rot_init, const Scalar epsilon) {
  // Total ops: 662

  // Unused inputs
  (void)encoder;
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (215)
  const Scalar _tmp0 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp1 = -2 * std::pow(_tmp0, Scalar(2));
  const Scalar _tmp2 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp3 = 1 - 2 * std::pow(_tmp2, Scalar(2));
  const Scalar _tmp4 = Scalar(0.20999999999999999) * _tmp1 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp5 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp6 = 2 * _tmp2;
  const Scalar _tmp7 = _tmp5 * _tmp6;
  const Scalar _tmp8 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp9 = _tmp0 * _tmp8;
  const Scalar _tmp10 =
      -Scalar(0.010999999999999999) * _tmp7 - Scalar(0.010999999999999999) * _tmp9;
  const Scalar _tmp11 = 2 * _tmp0 * _tmp5;
  const Scalar _tmp12 = _tmp2 * _tmp8;
  const Scalar _tmp13 = Scalar(0.20999999999999999) * _tmp11 - Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp14 = _tmp10 + _tmp13;
  const Scalar _tmp15 = _tmp14 + _tmp4;
  const Scalar _tmp16 = _tmp15 - p_c(0, 0) + position_vector(0, 0);
  const Scalar _tmp17 = -_tmp4;
  const Scalar _tmp18 = _tmp10 - _tmp13;
  const Scalar _tmp19 = _tmp17 + _tmp18;
  const Scalar _tmp20 = _tmp19 - p_a(0, 0) + position_vector(0, 0);
  const Scalar _tmp21 = Scalar(1.0) / (_tmp20);
  const Scalar _tmp22 = -2 * std::pow(_tmp5, Scalar(2));
  const Scalar _tmp23 = Scalar(0.20999999999999999) * _tmp22 + Scalar(0.20999999999999999) * _tmp3;
  const Scalar _tmp24 = -_tmp23;
  const Scalar _tmp25 = _tmp0 * _tmp6;
  const Scalar _tmp26 = _tmp5 * _tmp8;
  const Scalar _tmp27 =
      -Scalar(0.010999999999999999) * _tmp25 + Scalar(0.010999999999999999) * _tmp26;
  const Scalar _tmp28 = Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp12;
  const Scalar _tmp29 = _tmp27 - _tmp28;
  const Scalar _tmp30 = _tmp24 + _tmp29;
  const Scalar _tmp31 = _tmp30 - p_a(1, 0) + position_vector(1, 0);
  const Scalar _tmp32 = _tmp21 * _tmp31;
  const Scalar _tmp33 =
      std::sqrt(Scalar(std::pow(_tmp20, Scalar(2)) + std::pow(_tmp31, Scalar(2))));
  const Scalar _tmp34 = Scalar(1.0) / (_tmp33);
  const Scalar _tmp35 = _tmp21 * _tmp33;
  const Scalar _tmp36 = _tmp35 * (_tmp19 * _tmp31 * _tmp34 - _tmp20 * _tmp30 * _tmp34);
  const Scalar _tmp37 = -Scalar(0.010999999999999999) * _tmp1 -
                        Scalar(0.010999999999999999) * _tmp22 + Scalar(-0.010999999999999999);
  const Scalar _tmp38 = Scalar(0.20999999999999999) * _tmp7 - Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp39 = -_tmp38;
  const Scalar _tmp40 = Scalar(0.20999999999999999) * _tmp25 + Scalar(0.20999999999999999) * _tmp26;
  const Scalar _tmp41 = _tmp37 + _tmp39 - _tmp40;
  const Scalar _tmp42 = _tmp32 * _tmp41;
  const Scalar _tmp43 = _tmp37 + _tmp40;
  const Scalar _tmp44 = _tmp38 + _tmp43;
  const Scalar _tmp45 = _tmp27 + _tmp28;
  const Scalar _tmp46 = _tmp23 + _tmp45;
  const Scalar _tmp47 = _tmp46 - p_c(1, 0) + position_vector(1, 0);
  const Scalar _tmp48 = std::pow(_tmp16, Scalar(2));
  const Scalar _tmp49 = std::pow(_tmp47, Scalar(2));
  const Scalar _tmp50 = _tmp48 + _tmp49;
  const Scalar _tmp51 = std::pow(_tmp50, Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp52 = _tmp47 * _tmp51;
  const Scalar _tmp53 = _tmp41 * _tmp51;
  const Scalar _tmp54 = _tmp16 * _tmp53;
  const Scalar _tmp55 = -_tmp32 * _tmp54 + _tmp44 * _tmp52;
  const Scalar _tmp56 = _tmp32 * _tmp51;
  const Scalar _tmp57 = _tmp16 * _tmp56 - _tmp52;
  const Scalar _tmp58 = Scalar(1.0) / (_tmp57);
  const Scalar _tmp59 = _tmp32 * _tmp58;
  const Scalar _tmp60 = Scalar(1.0) * _tmp19;
  const Scalar _tmp61 = Scalar(1.0) * _tmp30;
  const Scalar _tmp62 = (-_tmp15 + _tmp60) / (_tmp46 - _tmp61);
  const Scalar _tmp63 = _tmp44 * _tmp51;
  const Scalar _tmp64 = -_tmp16 * _tmp63 + _tmp54;
  const Scalar _tmp65 = -_tmp41 + _tmp59 * _tmp64 - _tmp62 * (_tmp42 + _tmp55 * _tmp59);
  const Scalar _tmp66 = _tmp14 + _tmp17;
  const Scalar _tmp67 = _tmp66 - p_d(0, 0) + position_vector(0, 0);
  const Scalar _tmp68 = _tmp23 + _tmp29;
  const Scalar _tmp69 = _tmp68 - p_d(1, 0) + position_vector(1, 0);
  const Scalar _tmp70 = std::pow(Scalar(std::pow(_tmp67, Scalar(2)) + std::pow(_tmp69, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp71 = _tmp67 * _tmp70;
  const Scalar _tmp72 = _tmp41 * _tmp71;
  const Scalar _tmp73 = _tmp39 + _tmp43;
  const Scalar _tmp74 = _tmp69 * _tmp70;
  const Scalar _tmp75 = _tmp32 * _tmp71 - _tmp74;
  const Scalar _tmp76 = _tmp58 * _tmp75;
  const Scalar _tmp77 = -_tmp62 * (-_tmp32 * _tmp72 - _tmp55 * _tmp76 + _tmp73 * _tmp74) -
                        _tmp64 * _tmp76 - _tmp71 * _tmp73 + _tmp72;
  const Scalar _tmp78 = Scalar(1.0) / (_tmp77);
  const Scalar _tmp79 = _tmp46 * _tmp51;
  const Scalar _tmp80 = _tmp36 * _tmp51;
  const Scalar _tmp81 = -_tmp15 * _tmp52 + _tmp16 * _tmp79 + _tmp16 * _tmp80;
  const Scalar _tmp82 = _tmp36 * _tmp71 - _tmp66 * _tmp74 + _tmp68 * _tmp71 - _tmp76 * _tmp81;
  const Scalar _tmp83 = _tmp78 * _tmp82;
  const Scalar _tmp84 = -_tmp36 + _tmp59 * _tmp81 - _tmp65 * _tmp83;
  const Scalar _tmp85 = Scalar(1.0) / (_tmp82);
  const Scalar _tmp86 = _tmp77 * _tmp85;
  const Scalar _tmp87 = _tmp65 + _tmp84 * _tmp86;
  const Scalar _tmp88 = _tmp75 * _tmp78;
  const Scalar _tmp89 = -_tmp32 - _tmp87 * _tmp88;
  const Scalar _tmp90 = _tmp51 * _tmp58;
  const Scalar _tmp91 = _tmp89 * _tmp90;
  const Scalar _tmp92 = _tmp71 * _tmp78;
  const Scalar _tmp93 = _tmp18 + _tmp4;
  const Scalar _tmp94 = _tmp93 - p_b(0, 0) + position_vector(0, 0);
  const Scalar _tmp95 = _tmp24 + _tmp45;
  const Scalar _tmp96 = _tmp95 - p_b(1, 0) + position_vector(1, 0);
  const Scalar _tmp97 = std::pow(Scalar(std::pow(_tmp94, Scalar(2)) + std::pow(_tmp96, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp98 = _tmp94 * _tmp97;
  const Scalar _tmp99 = _tmp98 * fh1;
  const Scalar _tmp100 = _tmp35 * _tmp99;
  const Scalar _tmp101 = Scalar(1.0) * _tmp58;
  const Scalar _tmp102 = _tmp101 * _tmp62;
  const Scalar _tmp103 = -_tmp101 * _tmp64 + _tmp102 * _tmp55;
  const Scalar _tmp104 = -_tmp101 * _tmp81 - _tmp103 * _tmp83;
  const Scalar _tmp105 = _tmp103 + _tmp104 * _tmp86;
  const Scalar _tmp106 = -_tmp105 * _tmp88 + Scalar(1.0);
  const Scalar _tmp107 = _tmp106 * _tmp90;
  const Scalar _tmp108 = _tmp96 * _tmp97;
  const Scalar _tmp109 = _tmp108 * fh1;
  const Scalar _tmp110 = _tmp109 * _tmp35;
  const Scalar _tmp111 = _tmp60 + _tmp61 * _tmp62;
  const Scalar _tmp112 = 0;
  const Scalar _tmp113 = _tmp16 * _tmp51;
  const Scalar _tmp114 = _tmp112 * _tmp76;
  const Scalar _tmp115 = _tmp114 * _tmp78;
  const Scalar _tmp116 = Scalar(333.54000000000002) - fv1;
  const Scalar _tmp117 = _tmp116 * _tmp35;
  const Scalar _tmp118 = Scalar(1.0) * _tmp85;
  const Scalar _tmp119 = _tmp101 * _tmp75;
  const Scalar _tmp120 = _tmp119 * _tmp85;
  const Scalar _tmp121 = _tmp120 * _tmp51;
  const Scalar _tmp122 = fh1 * (_tmp108 * _tmp93 - _tmp95 * _tmp98);
  const Scalar _tmp123 = _tmp122 * _tmp35;
  const Scalar _tmp124 = std::exp(_tmp100 * (_tmp16 * _tmp91 + _tmp87 * _tmp92 + Scalar(1.0)) +
                                  _tmp110 * (_tmp105 * _tmp92 + _tmp107 * _tmp16) +
                                  _tmp117 * (_tmp112 * _tmp92 - _tmp113 * _tmp115) +
                                  _tmp123 * (_tmp118 * _tmp71 - _tmp121 * _tmp16));
  const Scalar _tmp125 = std::pow(_tmp50, Scalar(Scalar(-3) / Scalar(2)));
  const Scalar _tmp126 = _tmp125 * _tmp48;
  const Scalar _tmp127 = _tmp125 * _tmp16 * _tmp47;
  const Scalar _tmp128 = _tmp126 * _tmp36 + _tmp126 * _tmp46 - _tmp127 * _tmp15 - _tmp79 - _tmp80;
  const Scalar _tmp129 = _tmp126 * _tmp32 - _tmp127 - _tmp56;
  const Scalar _tmp130 = std::pow(_tmp57, Scalar(-2));
  const Scalar _tmp131 = _tmp130 * _tmp75;
  const Scalar _tmp132 = _tmp129 * _tmp131;
  const Scalar _tmp133 = -_tmp128 * _tmp76 + _tmp132 * _tmp81;
  const Scalar _tmp134 = std::pow(_tmp82, Scalar(-2));
  const Scalar _tmp135 = Scalar(1.0) * _tmp134;
  const Scalar _tmp136 = _tmp133 * _tmp135;
  const Scalar _tmp137 = _tmp113 * _tmp118;
  const Scalar _tmp138 = _tmp119 * _tmp134;
  const Scalar _tmp139 = _tmp113 * _tmp138;
  const Scalar _tmp140 = _tmp126 * _tmp58;
  const Scalar _tmp141 = _tmp130 * _tmp89;
  const Scalar _tmp142 = _tmp129 * _tmp141;
  const Scalar _tmp143 = _tmp126 * _tmp41;
  const Scalar _tmp144 = -_tmp126 * _tmp44 + _tmp143 - _tmp53 + _tmp63;
  const Scalar _tmp145 = _tmp127 * _tmp44;
  const Scalar _tmp146 = -_tmp143 * _tmp32 + _tmp145 + _tmp32 * _tmp53;
  const Scalar _tmp147 =
      _tmp132 * _tmp64 - _tmp144 * _tmp76 - _tmp62 * (_tmp132 * _tmp55 - _tmp146 * _tmp76);
  const Scalar _tmp148 = _tmp84 * _tmp85;
  const Scalar _tmp149 = _tmp134 * _tmp77;
  const Scalar _tmp150 = _tmp133 * _tmp149;
  const Scalar _tmp151 = _tmp130 * _tmp32;
  const Scalar _tmp152 = _tmp129 * _tmp151;
  const Scalar _tmp153 = _tmp65 * _tmp78;
  const Scalar _tmp154 =
      _tmp144 * _tmp59 - _tmp152 * _tmp64 - _tmp62 * (_tmp146 * _tmp59 - _tmp152 * _tmp55);
  const Scalar _tmp155 = std::pow(_tmp77, Scalar(-2));
  const Scalar _tmp156 = _tmp147 * _tmp155;
  const Scalar _tmp157 = _tmp156 * _tmp82;
  const Scalar _tmp158 = _tmp147 * _tmp148 - _tmp150 * _tmp84 + _tmp154 +
                         _tmp86 * (_tmp128 * _tmp59 - _tmp133 * _tmp153 - _tmp152 * _tmp81 -
                                   _tmp154 * _tmp83 + _tmp157 * _tmp65);
  const Scalar _tmp159 = _tmp75 * _tmp87;
  const Scalar _tmp160 = _tmp156 * _tmp159 - _tmp158 * _tmp88;
  const Scalar _tmp161 = _tmp16 * _tmp90;
  const Scalar _tmp162 = _tmp156 * _tmp71;
  const Scalar _tmp163 = _tmp112 * _tmp113 * _tmp130 * _tmp88;
  const Scalar _tmp164 = _tmp113 * _tmp114;
  const Scalar _tmp165 = _tmp103 * _tmp78;
  const Scalar _tmp166 = Scalar(1.0) * _tmp130;
  const Scalar _tmp167 = _tmp129 * _tmp166;
  const Scalar _tmp168 =
      -_tmp101 * _tmp144 + _tmp102 * _tmp146 - _tmp167 * _tmp55 * _tmp62 + _tmp167 * _tmp64;
  const Scalar _tmp169 = _tmp104 * _tmp85;
  const Scalar _tmp170 = -_tmp104 * _tmp150 + _tmp147 * _tmp169 + _tmp168 +
                         _tmp86 * (-_tmp101 * _tmp128 + _tmp103 * _tmp157 - _tmp133 * _tmp165 +
                                   _tmp167 * _tmp81 - _tmp168 * _tmp83);
  const Scalar _tmp171 = _tmp105 * _tmp71;
  const Scalar _tmp172 = _tmp106 * _tmp130;
  const Scalar _tmp173 = _tmp113 * _tmp172;
  const Scalar _tmp174 = _tmp105 * _tmp75;
  const Scalar _tmp175 = _tmp156 * _tmp174 - _tmp170 * _tmp88;
  const Scalar _tmp176 = _tmp106 * _tmp58;
  const Scalar _tmp177 = _tmp58 * _tmp99;
  const Scalar _tmp178 = _tmp112 * _tmp116;
  const Scalar _tmp179 = _tmp178 * _tmp78;
  const Scalar _tmp180 =
      std::exp(-_tmp109 * _tmp176 + _tmp120 * _tmp122 - _tmp177 * _tmp89 + _tmp179 * _tmp76);
  const Scalar _tmp181 = _tmp118 * _tmp122;
  const Scalar _tmp182 = _tmp156 * _tmp178;
  const Scalar _tmp183 = _tmp109 * _tmp58;
  const Scalar _tmp184 = _tmp122 * _tmp138;
  const Scalar _tmp185 = _tmp109 * _tmp172;
  const Scalar _tmp186 = _tmp78 * _tmp99;
  const Scalar _tmp187 = _tmp109 * _tmp78;
  const Scalar _tmp188 = std::exp(-_tmp105 * _tmp187 - _tmp179 - _tmp181 - _tmp186 * _tmp87);
  const Scalar _tmp189 = _tmp105 * _tmp109;
  const Scalar _tmp190 = _tmp87 * _tmp99;
  const Scalar _tmp191 = _tmp125 * _tmp49;
  const Scalar _tmp192 = _tmp127 * _tmp32 - _tmp191 + _tmp51;
  const Scalar _tmp193 = _tmp192 * _tmp64;
  const Scalar _tmp194 = -_tmp127 * _tmp42 + _tmp191 * _tmp44 - _tmp63;
  const Scalar _tmp195 = _tmp192 * _tmp55;
  const Scalar _tmp196 = _tmp127 * _tmp41 - _tmp145;
  const Scalar _tmp197 =
      _tmp131 * _tmp193 - _tmp196 * _tmp76 - _tmp62 * (_tmp131 * _tmp195 - _tmp194 * _tmp76);
  const Scalar _tmp198 = _tmp155 * _tmp197;
  const Scalar _tmp199 = _tmp198 * _tmp71;
  const Scalar _tmp200 = _tmp127 * _tmp36 + _tmp127 * _tmp46 - _tmp15 * _tmp191 + _tmp15 * _tmp51;
  const Scalar _tmp201 = _tmp192 * _tmp81;
  const Scalar _tmp202 = _tmp131 * _tmp201 - _tmp200 * _tmp76;
  const Scalar _tmp203 = _tmp149 * _tmp202;
  const Scalar _tmp204 = _tmp198 * _tmp82;
  const Scalar _tmp205 =
      -_tmp151 * _tmp193 + _tmp196 * _tmp59 - _tmp62 * (-_tmp151 * _tmp195 + _tmp194 * _tmp59);
  const Scalar _tmp206 = _tmp148 * _tmp197 - _tmp203 * _tmp84 + _tmp205 +
                         _tmp86 * (-_tmp151 * _tmp201 - _tmp153 * _tmp202 + _tmp200 * _tmp59 +
                                   _tmp204 * _tmp65 - _tmp205 * _tmp83);
  const Scalar _tmp207 = _tmp141 * _tmp192;
  const Scalar _tmp208 = _tmp159 * _tmp198 - _tmp206 * _tmp88;
  const Scalar _tmp209 =
      -_tmp101 * _tmp196 + _tmp102 * _tmp194 + _tmp166 * _tmp193 - _tmp166 * _tmp195 * _tmp62;
  const Scalar _tmp210 = -_tmp104 * _tmp203 + _tmp169 * _tmp197 + _tmp209 +
                         _tmp86 * (-_tmp101 * _tmp200 + _tmp103 * _tmp204 - _tmp165 * _tmp202 +
                                   _tmp166 * _tmp201 - _tmp209 * _tmp83);
  const Scalar _tmp211 = _tmp174 * _tmp198 - _tmp210 * _tmp88;
  const Scalar _tmp212 = _tmp135 * _tmp202;
  const Scalar _tmp213 = _tmp131 * _tmp192;
  const Scalar _tmp214 = _tmp178 * _tmp198;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 3> _res;

  _res(0, 0) = 0;
  _res(1, 0) = -_tmp124 * (-_tmp100 * (-_tmp113 * _tmp142 + _tmp140 * _tmp89 + _tmp158 * _tmp92 +
                                       _tmp160 * _tmp161 - _tmp162 * _tmp87 - _tmp91) -
                           _tmp110 * (_tmp106 * _tmp140 - _tmp107 - _tmp129 * _tmp173 -
                                      _tmp156 * _tmp171 + _tmp161 * _tmp175 + _tmp170 * _tmp92) -
                           _tmp117 * (-_tmp112 * _tmp162 - _tmp115 * _tmp126 + _tmp115 * _tmp51 +
                                      _tmp129 * _tmp163 + _tmp156 * _tmp164) -
                           _tmp123 * (-_tmp120 * _tmp126 + _tmp121 + _tmp132 * _tmp137 +
                                      _tmp133 * _tmp139 - _tmp136 * _tmp71));
  _res(2, 0) =
      -_tmp180 * (-_tmp129 * _tmp185 + _tmp132 * _tmp179 + _tmp132 * _tmp181 + _tmp133 * _tmp184 -
                  _tmp142 * _tmp99 + _tmp160 * _tmp177 + _tmp175 * _tmp183 + _tmp182 * _tmp76);
  _res(3, 0) = -_tmp188 * (-_tmp122 * _tmp136 - _tmp156 * _tmp189 - _tmp156 * _tmp190 +
                           _tmp158 * _tmp186 + _tmp170 * _tmp187 - _tmp182);
  _res(0, 1) = 0;
  _res(1, 1) =
      -_tmp124 *
      (-_tmp100 * (-_tmp113 * _tmp207 + _tmp127 * _tmp58 * _tmp89 + _tmp161 * _tmp208 -
                   _tmp199 * _tmp87 + _tmp206 * _tmp92) -
       _tmp110 * (_tmp127 * _tmp176 + _tmp161 * _tmp211 - _tmp171 * _tmp198 - _tmp173 * _tmp192 +
                  _tmp210 * _tmp92) -
       _tmp117 * (-_tmp112 * _tmp199 - _tmp115 * _tmp127 + _tmp163 * _tmp192 + _tmp164 * _tmp198) -
       _tmp123 * (-_tmp120 * _tmp127 + _tmp137 * _tmp213 + _tmp139 * _tmp202 - _tmp212 * _tmp71));
  _res(2, 1) =
      -_tmp180 * (_tmp177 * _tmp208 + _tmp179 * _tmp213 + _tmp181 * _tmp213 + _tmp183 * _tmp211 +
                  _tmp184 * _tmp202 - _tmp185 * _tmp192 - _tmp207 * _tmp99 + _tmp214 * _tmp76);
  _res(3, 1) = -_tmp188 * (-_tmp122 * _tmp212 + _tmp186 * _tmp206 + _tmp187 * _tmp210 -
                           _tmp189 * _tmp198 - _tmp190 * _tmp198 - _tmp214);
  _res(0, 2) = 0;
  _res(1, 2) = 0;
  _res(2, 2) = 0;
  _res(3, 2) = 0;

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
