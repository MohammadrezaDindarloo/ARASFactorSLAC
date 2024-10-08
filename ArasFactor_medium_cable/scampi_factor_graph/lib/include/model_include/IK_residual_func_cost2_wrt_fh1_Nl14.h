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
 * Symbolic function: IK_residual_func_cost2_wrt_fh1_Nl14
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost2WrtFh1Nl14(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 646

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (216)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp2 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp3 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp4 = 1 - 2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp4;
  const Scalar _tmp6 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp7 = 2 * _tmp1 * _tmp6;
  const Scalar _tmp8 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                       2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp9 = _tmp3 * _tmp8;
  const Scalar _tmp10 = _tmp7 + _tmp9;
  const Scalar _tmp11 = -Scalar(0.010999999999999999) * _tmp10;
  const Scalar _tmp12 = 2 * _tmp3;
  const Scalar _tmp13 = _tmp12 * _tmp6;
  const Scalar _tmp14 = _tmp1 * _tmp8;
  const Scalar _tmp15 = Scalar(0.20999999999999999) * _tmp13 - Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp16 = _tmp11 + _tmp15;
  const Scalar _tmp17 = _tmp16 + _tmp5;
  const Scalar _tmp18 = _tmp17 + position_vector(0, 0);
  const Scalar _tmp19 = Scalar(0.20999999999999999) * _tmp13 + Scalar(0.20999999999999999) * _tmp14;
  const Scalar _tmp20 = -2 * std::pow(_tmp6, Scalar(2));
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp20 +
                        Scalar(0.20999999999999999);
  const Scalar _tmp22 = _tmp1 * _tmp12;
  const Scalar _tmp23 = _tmp6 * _tmp8;
  const Scalar _tmp24 = _tmp22 - _tmp23;
  const Scalar _tmp25 = -Scalar(0.010999999999999999) * _tmp24;
  const Scalar _tmp26 = _tmp21 + _tmp25;
  const Scalar _tmp27 = _tmp19 + _tmp26;
  const Scalar _tmp28 = _tmp27 + position_vector(1, 0);
  const Scalar _tmp29 = _tmp0 * fv1;
  const Scalar _tmp30 = std::asinh(_tmp29);
  const Scalar _tmp31 = Scalar(9.6622558468725703) * _tmp30;
  const Scalar _tmp32 =
      -Scalar(0.1034955) * _tmp31 * fh1 -
      Scalar(0.49413274378274363) *
          std::sqrt(
              Scalar(Scalar(0.32387954179207445) *
                         std::pow(Scalar(1 - Scalar(0.36803241838814449) * _tmp18), Scalar(2)) +
                     std::pow(Scalar(1 - Scalar(0.2094487793051498) * _tmp28), Scalar(2))));
  const Scalar _tmp33 = _tmp0 * _tmp32;
  const Scalar _tmp34 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp35 =
      std::pow(Scalar(_tmp34 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp36 = Scalar(1.0) * _tmp30;
  const Scalar _tmp37 = -_tmp21 + _tmp25;
  const Scalar _tmp38 = _tmp19 + _tmp37;
  const Scalar _tmp39 = _tmp38 + position_vector(1, 0);
  const Scalar _tmp40 = _tmp11 - _tmp15;
  const Scalar _tmp41 = _tmp40 + _tmp5;
  const Scalar _tmp42 = _tmp41 + position_vector(0, 0);
  const Scalar _tmp43 = -_tmp5;
  const Scalar _tmp44 = _tmp16 + _tmp43;
  const Scalar _tmp45 = _tmp44 + position_vector(0, 0);
  const Scalar _tmp46 = _tmp45 + Scalar(1.7965602546229);
  const Scalar _tmp47 = -_tmp19;
  const Scalar _tmp48 = _tmp26 + _tmp47;
  const Scalar _tmp49 = _tmp48 + position_vector(1, 0);
  const Scalar _tmp50 = _tmp49 + Scalar(-4.83288938413423);
  const Scalar _tmp51 = std::pow(Scalar(std::pow(_tmp46, Scalar(2)) + std::pow(_tmp50, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp52 = _tmp46 * _tmp51;
  const Scalar _tmp53 = Scalar(0.20999999999999999) * _tmp22 + Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp54 = -_tmp53;
  const Scalar _tmp55 =
      -Scalar(0.010999999999999999) * _tmp20 - Scalar(0.010999999999999999) * _tmp4;
  const Scalar _tmp56 = Scalar(0.20999999999999999) * _tmp7 - Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp57 = _tmp55 + _tmp56;
  const Scalar _tmp58 = _tmp54 + _tmp57;
  const Scalar _tmp59 = _tmp42 + Scalar(-2.5193355532036801);
  const Scalar _tmp60 = Scalar(1.0) / (_tmp59);
  const Scalar _tmp61 = _tmp39 + Scalar(8.3885017487099702);
  const Scalar _tmp62 = _tmp60 * _tmp61;
  const Scalar _tmp63 = _tmp58 * _tmp62;
  const Scalar _tmp64 = _tmp55 - _tmp56;
  const Scalar _tmp65 = _tmp54 + _tmp64;
  const Scalar _tmp66 = _tmp37 + _tmp47;
  const Scalar _tmp67 = _tmp66 + position_vector(1, 0);
  const Scalar _tmp68 = _tmp67 + Scalar(8.3196563720703107);
  const Scalar _tmp69 = _tmp40 + _tmp43;
  const Scalar _tmp70 = _tmp69 + position_vector(0, 0);
  const Scalar _tmp71 = _tmp70 + Scalar(1.9874742031097401);
  const Scalar _tmp72 = std::pow(Scalar(std::pow(_tmp68, Scalar(2)) + std::pow(_tmp71, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp73 = _tmp68 * _tmp72;
  const Scalar _tmp74 = _tmp71 * _tmp72;
  const Scalar _tmp75 = -_tmp63 * _tmp74 + _tmp65 * _tmp73;
  const Scalar _tmp76 = _tmp50 * _tmp51;
  const Scalar _tmp77 = _tmp52 * _tmp62 - _tmp76;
  const Scalar _tmp78 = Scalar(1.0) / (_tmp62 * _tmp74 - _tmp73);
  const Scalar _tmp79 = _tmp77 * _tmp78;
  const Scalar _tmp80 = _tmp53 + _tmp64;
  const Scalar _tmp81 = -_tmp52 * _tmp63 - _tmp75 * _tmp79 + _tmp76 * _tmp80;
  const Scalar _tmp82 = Scalar(1.0) * _tmp38;
  const Scalar _tmp83 = -_tmp82;
  const Scalar _tmp84 = Scalar(1.0) / (_tmp66 + _tmp83);
  const Scalar _tmp85 = Scalar(1.0) * _tmp41;
  const Scalar _tmp86 = -_tmp69 + _tmp85;
  const Scalar _tmp87 = _tmp84 * _tmp86;
  const Scalar _tmp88 = _tmp78 * (_tmp58 * _tmp74 - _tmp65 * _tmp74);
  const Scalar _tmp89 = _tmp52 * _tmp58 - _tmp52 * _tmp80 - _tmp77 * _tmp88 - _tmp81 * _tmp87;
  const Scalar _tmp90 = Scalar(1.0) / (_tmp89);
  const Scalar _tmp91 =
      std::sqrt(Scalar(std::pow(_tmp59, Scalar(2)) + std::pow(_tmp61, Scalar(2))));
  const Scalar _tmp92 = Scalar(1.0) / (_tmp91);
  const Scalar _tmp93 = _tmp60 * _tmp91;
  const Scalar _tmp94 = _tmp93 * (-_tmp38 * _tmp59 * _tmp92 + _tmp41 * _tmp61 * _tmp92);
  const Scalar _tmp95 = _tmp66 * _tmp74 - _tmp69 * _tmp73 + _tmp74 * _tmp94;
  const Scalar _tmp96 = _tmp62 * _tmp78;
  const Scalar _tmp97 = _tmp63 + _tmp75 * _tmp96;
  const Scalar _tmp98 = -_tmp58 + _tmp62 * _tmp88 - _tmp87 * _tmp97;
  const Scalar _tmp99 = -_tmp44 * _tmp76 + _tmp48 * _tmp52 + _tmp52 * _tmp94 - _tmp79 * _tmp95;
  const Scalar _tmp100 = _tmp90 * _tmp99;
  const Scalar _tmp101 = Scalar(1.0) / (_tmp99);
  const Scalar _tmp102 = _tmp101 * _tmp89;
  const Scalar _tmp103 = _tmp102 * (-_tmp100 * _tmp98 - _tmp94 + _tmp95 * _tmp96);
  const Scalar _tmp104 = _tmp103 + _tmp98;
  const Scalar _tmp105 = _tmp104 * _tmp90;
  const Scalar _tmp106 = _tmp77 * _tmp90;
  const Scalar _tmp107 = _tmp78 * (-_tmp104 * _tmp106 - _tmp62);
  const Scalar _tmp108 = _tmp18 + Scalar(-2.7171519410699099);
  const Scalar _tmp109 = _tmp28 + Scalar(-4.7744369927459998);
  const Scalar _tmp110 =
      std::pow(Scalar(std::pow(_tmp108, Scalar(2)) + std::pow(_tmp109, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp111 = _tmp108 * _tmp110;
  const Scalar _tmp112 = _tmp111 * _tmp93 * (_tmp105 * _tmp52 + _tmp107 * _tmp74 + Scalar(1.0));
  const Scalar _tmp113 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp114 = _tmp82 * _tmp87 + _tmp85;
  const Scalar _tmp115 = 0;
  const Scalar _tmp116 = _tmp115 * _tmp90;
  const Scalar _tmp117 = _tmp74 * _tmp79;
  const Scalar _tmp118 = _tmp109 * _tmp110;
  const Scalar _tmp119 = -_tmp111 * _tmp27 + _tmp118 * _tmp17;
  const Scalar _tmp120 = Scalar(1.0) * _tmp101;
  const Scalar _tmp121 = _tmp119 * _tmp93 * (-_tmp117 * _tmp120 + _tmp120 * _tmp52);
  const Scalar _tmp122 = Scalar(1.0) * _tmp84;
  const Scalar _tmp123 = _tmp122 * _tmp75 * _tmp78 * _tmp86 - Scalar(1.0) * _tmp88;
  const Scalar _tmp124 = Scalar(1.0) * _tmp78;
  const Scalar _tmp125 = _tmp102 * (-_tmp100 * _tmp123 - _tmp124 * _tmp95);
  const Scalar _tmp126 = _tmp123 + _tmp125;
  const Scalar _tmp127 = _tmp78 * (-_tmp106 * _tmp126 + Scalar(1.0));
  const Scalar _tmp128 = _tmp126 * _tmp90;
  const Scalar _tmp129 = _tmp118 * _tmp93 * (_tmp127 * _tmp74 + _tmp128 * _tmp52);
  const Scalar _tmp130 = -_tmp112 * fh1 -
                         _tmp113 * _tmp93 * (-_tmp116 * _tmp117 + _tmp116 * _tmp52) -
                         _tmp121 * fh1 - _tmp129 * fh1;
  const Scalar _tmp131 = Scalar(1.0) / (_tmp130);
  const Scalar _tmp132 = _tmp53 + _tmp57;
  const Scalar _tmp133 = _tmp132 * fh1;
  const Scalar _tmp134 = -_tmp118 * _tmp133 - Scalar(5.1796800000000003) * _tmp24 - _tmp27 * fv1;
  const Scalar _tmp135 = _tmp48 + _tmp83;
  const Scalar _tmp136 = _tmp135 * _tmp87;
  const Scalar _tmp137 = Scalar(1.0) / (-_tmp136 - _tmp44 + _tmp85);
  const Scalar _tmp138 = Scalar(1.0) * _tmp137;
  const Scalar _tmp139 = _tmp136 * _tmp138 + Scalar(1.0);
  const Scalar _tmp140 = _tmp138 * _tmp87;
  const Scalar _tmp141 = -Scalar(1.0) * _tmp122 * _tmp139 + Scalar(1.0) * _tmp140;
  const Scalar _tmp142 = _tmp135 * _tmp137;
  const Scalar _tmp143 = _tmp81 * _tmp90;
  const Scalar _tmp144 = _tmp103 * _tmp142 - _tmp104 * _tmp143 + _tmp97;
  const Scalar _tmp145 = Scalar(1.0) * _tmp111 * (_tmp103 * _tmp138 - _tmp122 * _tmp144);
  const Scalar _tmp146 = _tmp114 * _tmp137;
  const Scalar _tmp147 = _tmp84 * (-_tmp115 * _tmp143 - _tmp135 * _tmp146 + _tmp83);
  const Scalar _tmp148 = -_tmp124 * _tmp75 + _tmp125 * _tmp142 - _tmp126 * _tmp143;
  const Scalar _tmp149 = Scalar(1.0) * _tmp118 * (-_tmp122 * _tmp148 + _tmp125 * _tmp138);
  const Scalar _tmp150 = _tmp102 * _tmp138;
  const Scalar _tmp151 = _tmp84 * (-_tmp120 * _tmp81 + _tmp135 * _tmp150);
  const Scalar _tmp152 = Scalar(1.0) * _tmp119;
  const Scalar _tmp153 = _tmp152 * (_tmp150 - Scalar(1.0) * _tmp151);
  const Scalar _tmp154 = Scalar(5.1796800000000003) * _tmp10 + _tmp111 * _tmp133 + _tmp17 * fv1;
  const Scalar _tmp155 = _tmp135 * _tmp84;
  const Scalar _tmp156 = Scalar(1.0) * _tmp138 * _tmp155 - Scalar(1.0) * _tmp138;
  const Scalar _tmp157 =
      Scalar(1.0) * _tmp113 * (-_tmp114 * _tmp138 - Scalar(1.0) * _tmp147 + Scalar(1.0)) +
      _tmp134 * _tmp141 + _tmp145 * fh1 + _tmp149 * fh1 + _tmp153 * fh1 + _tmp154 * _tmp156;
  const Scalar _tmp158 = std::asinh(_tmp131 * _tmp157);
  const Scalar _tmp159 = Scalar(9.6622558468725703) * _tmp130;
  const Scalar _tmp160 =
      -_tmp158 * _tmp159 -
      Scalar(8.3885017487099702) *
          std::sqrt(
              Scalar(Scalar(0.090199313518583735) *
                         std::pow(Scalar(1 - Scalar(0.39693005512043167) * _tmp42), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.11921079949155229) * _tmp39 - 1), Scalar(2))));
  const Scalar _tmp161 = Scalar(0.1034955) * _tmp131;
  const Scalar _tmp162 = _tmp160 * _tmp161;
  const Scalar _tmp163 = std::pow(_tmp130, Scalar(-2));
  const Scalar _tmp164 = _tmp111 * _tmp132;
  const Scalar _tmp165 = _tmp118 * _tmp132;
  const Scalar _tmp166 = -_tmp112 - _tmp121 - _tmp129;
  const Scalar _tmp167 = _tmp163 * _tmp166;
  const Scalar _tmp168 =
      (_tmp131 * (-_tmp141 * _tmp165 + _tmp145 + _tmp149 + _tmp153 + _tmp156 * _tmp164) -
       _tmp157 * _tmp167) /
      std::sqrt(Scalar(std::pow(_tmp157, Scalar(2)) * _tmp163 + 1));
  const Scalar _tmp169 = Scalar(9.6622558468725703) * _tmp166;
  const Scalar _tmp170 = Scalar(1.0) * _tmp158;
  const Scalar _tmp171 = _tmp113 * _tmp116;
  const Scalar _tmp172 = _tmp107 * _tmp111;
  const Scalar _tmp173 = _tmp118 * _tmp127;
  const Scalar _tmp174 = _tmp101 * _tmp152;
  const Scalar _tmp175 = _tmp174 * fh1;
  const Scalar _tmp176 = -_tmp171 * _tmp79 + _tmp172 * fh1 + _tmp173 * fh1 - _tmp175 * _tmp79;
  const Scalar _tmp177 = Scalar(1.0) / (_tmp176);
  const Scalar _tmp178 = _tmp134 * _tmp84;
  const Scalar _tmp179 = _tmp118 * _tmp84;
  const Scalar _tmp180 = _tmp148 * _tmp179;
  const Scalar _tmp181 = _tmp119 * _tmp151;
  const Scalar _tmp182 = _tmp111 * _tmp84;
  const Scalar _tmp183 = _tmp144 * _tmp182;
  const Scalar _tmp184 = _tmp138 * _tmp154;
  const Scalar _tmp185 = _tmp113 * _tmp147 + _tmp139 * _tmp178 - _tmp155 * _tmp184 + _tmp180 * fh1 +
                         _tmp181 * fh1 + _tmp183 * fh1;
  const Scalar _tmp186 = std::asinh(_tmp177 * _tmp185);
  const Scalar _tmp187 = Scalar(1.0) * _tmp186;
  const Scalar _tmp188 = std::pow(_tmp176, Scalar(-2));
  const Scalar _tmp189 = _tmp172 + _tmp173 - _tmp174 * _tmp79;
  const Scalar _tmp190 = _tmp188 * _tmp189;
  const Scalar _tmp191 = (_tmp177 * (-_tmp132 * _tmp135 * _tmp138 * _tmp182 -
                                     _tmp132 * _tmp139 * _tmp179 + _tmp180 + _tmp181 + _tmp183) -
                          _tmp185 * _tmp190) /
                         std::sqrt(Scalar(std::pow(_tmp185, Scalar(2)) * _tmp188 + 1));
  const Scalar _tmp192 = Scalar(9.6622558468725703) * _tmp176;
  const Scalar _tmp193 = Scalar(9.6622558468725703) * _tmp189;
  const Scalar _tmp194 = Scalar(0.1034955) * _tmp177;
  const Scalar _tmp195 =
      -_tmp186 * _tmp192 -
      Scalar(8.3196563720703107) *
          std::sqrt(
              Scalar(std::pow(Scalar(-Scalar(0.12019727201198747) * _tmp67 - 1), Scalar(2)) +
                     Scalar(0.057067943527034905) *
                         std::pow(Scalar(-Scalar(0.50315118477277876) * _tmp70 - 1), Scalar(2))));
  const Scalar _tmp196 = _tmp194 * _tmp195;
  const Scalar _tmp197 = _tmp118 * _tmp125 * _tmp137;
  const Scalar _tmp198 = _tmp103 * _tmp111 * _tmp137;
  const Scalar _tmp199 = _tmp119 * _tmp150;
  const Scalar _tmp200 = _tmp113 * _tmp146 - _tmp138 * _tmp178 * _tmp86 + _tmp184 - _tmp197 * fh1 -
                         _tmp198 * fh1 - _tmp199 * fh1;
  const Scalar _tmp201 = _tmp105 * _tmp111;
  const Scalar _tmp202 = _tmp118 * _tmp128;
  const Scalar _tmp203 = _tmp171 + _tmp175 + _tmp201 * fh1 + _tmp202 * fh1;
  const Scalar _tmp204 = Scalar(1.0) / (_tmp203);
  const Scalar _tmp205 = std::asinh(_tmp200 * _tmp204);
  const Scalar _tmp206 = Scalar(1.0) * _tmp205;
  const Scalar _tmp207 = std::pow(_tmp203, Scalar(-2));
  const Scalar _tmp208 = _tmp174 + _tmp201 + _tmp202;
  const Scalar _tmp209 = _tmp207 * _tmp208;
  const Scalar _tmp210 = (-_tmp200 * _tmp209 + _tmp204 * (_tmp138 * _tmp164 + _tmp140 * _tmp165 -
                                                          _tmp197 - _tmp198 - _tmp199)) /
                         std::sqrt(Scalar(std::pow(_tmp200, Scalar(2)) * _tmp207 + 1));
  const Scalar _tmp211 = Scalar(9.6622558468725703) * _tmp203;
  const Scalar _tmp212 =
      -_tmp205 * _tmp211 -
      Scalar(4.83288938413423) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.20691555724053495) * _tmp49), Scalar(2)) +
                     Scalar(0.13818785160942856) *
                         std::pow(Scalar(-Scalar(0.55661923802822921) * _tmp45 - 1), Scalar(2))));
  const Scalar _tmp213 = Scalar(0.1034955) * _tmp204;
  const Scalar _tmp214 = _tmp212 * _tmp213;
  const Scalar _tmp215 = Scalar(9.6622558468725703) * _tmp208;

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) =
      Scalar(9.6622558468725703) * fh1 *
          (Scalar(1.0) * _tmp34 * _tmp35 * fv1 * std::cosh(_tmp36) -
           (Scalar(0.1034955) * _tmp0 * (Scalar(9.6622558468725703) * _tmp29 * _tmp35 - _tmp31) -
            _tmp32 * _tmp34) *
               std::cosh(_tmp33)) -
      Scalar(9.6622558468725703) * std::sinh(_tmp33) -
      Scalar(9.6622558468725703) * std::sinh(_tmp36);
  _res(1, 0) = _tmp159 * (-Scalar(1.0) * _tmp168 * std::cosh(_tmp170) -
                          (-Scalar(0.1034955) * _tmp160 * _tmp167 +
                           _tmp161 * (-_tmp158 * _tmp169 - _tmp159 * _tmp168)) *
                              std::cosh(_tmp162)) +
               _tmp169 * (-std::sinh(_tmp162) - std::sinh(_tmp170));
  _res(2, 0) = _tmp192 * (-Scalar(1.0) * _tmp191 * std::cosh(_tmp187) -
                          (-Scalar(0.1034955) * _tmp190 * _tmp195 +
                           _tmp194 * (-_tmp186 * _tmp193 - _tmp191 * _tmp192)) *
                              std::cosh(_tmp196)) +
               _tmp193 * (-std::sinh(_tmp187) - std::sinh(_tmp196));
  _res(3, 0) = _tmp211 * (-Scalar(1.0) * _tmp210 * std::cosh(_tmp206) -
                          (-Scalar(0.1034955) * _tmp209 * _tmp212 +
                           _tmp213 * (-_tmp205 * _tmp215 - _tmp210 * _tmp211)) *
                              std::cosh(_tmp214)) +
               _tmp215 * (-std::sinh(_tmp206) - std::sinh(_tmp214));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
