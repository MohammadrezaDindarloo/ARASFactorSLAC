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
 * Symbolic function: IK_residual_func_cost1_wrt_fh1_Nl9
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
Eigen::Matrix<Scalar, 4, 1> IkResidualFuncCost1WrtFh1Nl9(
    const Scalar fh1, const Scalar fv1, const sym::Rot3<Scalar>& DeltaRot,
    const Eigen::Matrix<Scalar, 3, 1>& position_vector, const sym::Rot3<Scalar>& Rot_init,
    const Scalar epsilon) {
  // Total ops: 653

  // Unused inputs
  (void)epsilon;

  // Input arrays
  const Eigen::Matrix<Scalar, 4, 1>& _DeltaRot = DeltaRot.Data();
  const Eigen::Matrix<Scalar, 4, 1>& _Rot_init = Rot_init.Data();

  // Intermediate terms (210)
  const Scalar _tmp0 = Scalar(1.0) / (fh1);
  const Scalar _tmp1 = _DeltaRot[0] * _Rot_init[3] - _DeltaRot[1] * _Rot_init[2] +
                       _DeltaRot[2] * _Rot_init[1] + _DeltaRot[3] * _Rot_init[0];
  const Scalar _tmp2 = -2 * std::pow(_tmp1, Scalar(2));
  const Scalar _tmp3 = -_DeltaRot[0] * _Rot_init[1] + _DeltaRot[1] * _Rot_init[0] +
                       _DeltaRot[2] * _Rot_init[3] + _DeltaRot[3] * _Rot_init[2];
  const Scalar _tmp4 = -2 * std::pow(_tmp3, Scalar(2));
  const Scalar _tmp5 = Scalar(0.20999999999999999) * _tmp2 + Scalar(0.20999999999999999) * _tmp4 +
                       Scalar(0.20999999999999999);
  const Scalar _tmp6 = -_tmp5;
  const Scalar _tmp7 = _DeltaRot[0] * _Rot_init[2] + _DeltaRot[1] * _Rot_init[3] -
                       _DeltaRot[2] * _Rot_init[0] + _DeltaRot[3] * _Rot_init[1];
  const Scalar _tmp8 = 2 * _tmp3;
  const Scalar _tmp9 = _tmp7 * _tmp8;
  const Scalar _tmp10 = -2 * _DeltaRot[0] * _Rot_init[0] - 2 * _DeltaRot[1] * _Rot_init[1] -
                        2 * _DeltaRot[2] * _Rot_init[2] + 2 * _DeltaRot[3] * _Rot_init[3];
  const Scalar _tmp11 = _tmp1 * _tmp10;
  const Scalar _tmp12 = -_tmp11 + _tmp9;
  const Scalar _tmp13 = -Scalar(0.010999999999999999) * _tmp12;
  const Scalar _tmp14 = 2 * _tmp1 * _tmp7;
  const Scalar _tmp15 = _tmp10 * _tmp3;
  const Scalar _tmp16 = Scalar(0.20999999999999999) * _tmp14 + Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp17 = _tmp13 + _tmp16;
  const Scalar _tmp18 = _tmp17 + _tmp6;
  const Scalar _tmp19 = _tmp18 + position_vector(1, 0);
  const Scalar _tmp20 = 1 - 2 * std::pow(_tmp7, Scalar(2));
  const Scalar _tmp21 = Scalar(0.20999999999999999) * _tmp20 + Scalar(0.20999999999999999) * _tmp4;
  const Scalar _tmp22 = _tmp1 * _tmp8;
  const Scalar _tmp23 = _tmp10 * _tmp7;
  const Scalar _tmp24 = _tmp22 + _tmp23;
  const Scalar _tmp25 = -Scalar(0.010999999999999999) * _tmp24;
  const Scalar _tmp26 = Scalar(0.20999999999999999) * _tmp14 - Scalar(0.20999999999999999) * _tmp15;
  const Scalar _tmp27 = _tmp25 - _tmp26;
  const Scalar _tmp28 = _tmp21 + _tmp27;
  const Scalar _tmp29 = _tmp28 + position_vector(0, 0);
  const Scalar _tmp30 = _tmp0 * fv1;
  const Scalar _tmp31 = std::asinh(_tmp30);
  const Scalar _tmp32 = Scalar(9.6622558468725703) * _tmp31;
  const Scalar _tmp33 =
      -_tmp32 * fh1 -
      Scalar(8.3885017487099702) *
          std::sqrt(
              Scalar(Scalar(0.090199313518583735) *
                         std::pow(Scalar(1 - Scalar(0.39693005512043167) * _tmp29), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.11921079949155229) * _tmp19 - 1), Scalar(2))));
  const Scalar _tmp34 = Scalar(0.1034955) * _tmp0;
  const Scalar _tmp35 = _tmp33 * _tmp34;
  const Scalar _tmp36 = std::pow(fh1, Scalar(-2));
  const Scalar _tmp37 =
      std::pow(Scalar(_tmp36 * std::pow(fv1, Scalar(2)) + 1), Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp38 = Scalar(1.0) * _tmp31;
  const Scalar _tmp39 = -_tmp21;
  const Scalar _tmp40 = _tmp27 + _tmp39;
  const Scalar _tmp41 = _tmp40 + position_vector(0, 0);
  const Scalar _tmp42 = _tmp41 + Scalar(1.9874742031097401);
  const Scalar _tmp43 = _tmp13 - _tmp16;
  const Scalar _tmp44 = _tmp43 + _tmp6;
  const Scalar _tmp45 = _tmp44 + position_vector(1, 0);
  const Scalar _tmp46 = _tmp45 + Scalar(8.3196563720703107);
  const Scalar _tmp47 = std::pow(Scalar(std::pow(_tmp42, Scalar(2)) + std::pow(_tmp46, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp48 = _tmp42 * _tmp47;
  const Scalar _tmp49 = Scalar(0.20999999999999999) * _tmp22 - Scalar(0.20999999999999999) * _tmp23;
  const Scalar _tmp50 = -_tmp49;
  const Scalar _tmp51 =
      -Scalar(0.010999999999999999) * _tmp2 - Scalar(0.010999999999999999) * _tmp20;
  const Scalar _tmp52 = Scalar(0.20999999999999999) * _tmp11 + Scalar(0.20999999999999999) * _tmp9;
  const Scalar _tmp53 = _tmp51 - _tmp52;
  const Scalar _tmp54 = _tmp50 + _tmp53;
  const Scalar _tmp55 = _tmp51 + _tmp52;
  const Scalar _tmp56 = _tmp49 + _tmp55;
  const Scalar _tmp57 = _tmp25 + _tmp26;
  const Scalar _tmp58 = _tmp39 + _tmp57;
  const Scalar _tmp59 = _tmp58 + position_vector(0, 0);
  const Scalar _tmp60 = _tmp59 + Scalar(1.7965602546229);
  const Scalar _tmp61 = _tmp43 + _tmp5;
  const Scalar _tmp62 = _tmp61 + position_vector(1, 0);
  const Scalar _tmp63 = _tmp62 + Scalar(-4.83288938413423);
  const Scalar _tmp64 = std::pow(Scalar(std::pow(_tmp60, Scalar(2)) + std::pow(_tmp63, Scalar(2))),
                                 Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp65 = _tmp60 * _tmp64;
  const Scalar _tmp66 = _tmp50 + _tmp55;
  const Scalar _tmp67 = _tmp56 * _tmp65 - _tmp65 * _tmp66;
  const Scalar _tmp68 = _tmp46 * _tmp47;
  const Scalar _tmp69 = _tmp17 + _tmp5;
  const Scalar _tmp70 = _tmp69 + position_vector(1, 0);
  const Scalar _tmp71 = _tmp70 + Scalar(-4.7744369927459998);
  const Scalar _tmp72 = _tmp21 + _tmp57;
  const Scalar _tmp73 = _tmp72 + position_vector(0, 0);
  const Scalar _tmp74 = _tmp73 + Scalar(-2.7171519410699099);
  const Scalar _tmp75 = Scalar(1.0) / (_tmp74);
  const Scalar _tmp76 = _tmp71 * _tmp75;
  const Scalar _tmp77 = _tmp48 * _tmp76 - _tmp68;
  const Scalar _tmp78 = _tmp63 * _tmp64;
  const Scalar _tmp79 = Scalar(1.0) / (_tmp65 * _tmp76 - _tmp78);
  const Scalar _tmp80 = _tmp77 * _tmp79;
  const Scalar _tmp81 = _tmp56 * _tmp76;
  const Scalar _tmp82 = -_tmp65 * _tmp81 + _tmp66 * _tmp78;
  const Scalar _tmp83 = -_tmp48 * _tmp81 + _tmp54 * _tmp68 - _tmp80 * _tmp82;
  const Scalar _tmp84 = Scalar(1.0) * _tmp69;
  const Scalar _tmp85 = -_tmp84;
  const Scalar _tmp86 = Scalar(1.0) / (_tmp61 + _tmp85);
  const Scalar _tmp87 = Scalar(1.0) * _tmp72;
  const Scalar _tmp88 = _tmp86 * (-_tmp58 + _tmp87);
  const Scalar _tmp89 = -_tmp48 * _tmp54 + _tmp48 * _tmp56 - _tmp67 * _tmp80 - _tmp83 * _tmp88;
  const Scalar _tmp90 = Scalar(1.0) / (_tmp89);
  const Scalar _tmp91 = _tmp76 * _tmp79;
  const Scalar _tmp92 = _tmp81 + _tmp82 * _tmp91;
  const Scalar _tmp93 = -_tmp56 + _tmp67 * _tmp91 - _tmp88 * _tmp92;
  const Scalar _tmp94 =
      std::sqrt(Scalar(std::pow(_tmp71, Scalar(2)) + std::pow(_tmp74, Scalar(2))));
  const Scalar _tmp95 = Scalar(1.0) / (_tmp94);
  const Scalar _tmp96 = _tmp75 * _tmp94;
  const Scalar _tmp97 = _tmp96 * (-_tmp69 * _tmp74 * _tmp95 + _tmp71 * _tmp72 * _tmp95);
  const Scalar _tmp98 = _tmp79 * (-_tmp58 * _tmp78 + _tmp61 * _tmp65 + _tmp65 * _tmp97);
  const Scalar _tmp99 = -_tmp40 * _tmp68 + _tmp44 * _tmp48 + _tmp48 * _tmp97 - _tmp77 * _tmp98;
  const Scalar _tmp100 = _tmp90 * _tmp99;
  const Scalar _tmp101 = Scalar(1.0) / (_tmp99);
  const Scalar _tmp102 = _tmp101 * _tmp89;
  const Scalar _tmp103 = _tmp102 * (-_tmp100 * _tmp93 + _tmp76 * _tmp98 - _tmp97);
  const Scalar _tmp104 = _tmp90 * (_tmp103 + _tmp93);
  const Scalar _tmp105 = _tmp79 * (-_tmp104 * _tmp77 - _tmp76);
  const Scalar _tmp106 = _tmp29 + Scalar(-2.5193355532036801);
  const Scalar _tmp107 = _tmp19 + Scalar(8.3885017487099702);
  const Scalar _tmp108 =
      std::pow(Scalar(std::pow(_tmp106, Scalar(2)) + std::pow(_tmp107, Scalar(2))),
               Scalar(Scalar(-1) / Scalar(2)));
  const Scalar _tmp109 = _tmp106 * _tmp108;
  const Scalar _tmp110 = _tmp109 * _tmp96 * (_tmp104 * _tmp48 + _tmp105 * _tmp65 + Scalar(1.0));
  const Scalar _tmp111 = Scalar(1.0) * _tmp79;
  const Scalar _tmp112 = _tmp111 * _tmp82;
  const Scalar _tmp113 = -_tmp111 * _tmp67 + _tmp112 * _tmp88;
  const Scalar _tmp114 = _tmp102 * (-_tmp100 * _tmp113 - Scalar(1.0) * _tmp98);
  const Scalar _tmp115 = _tmp90 * (_tmp113 + _tmp114);
  const Scalar _tmp116 = _tmp79 * (-_tmp115 * _tmp77 + Scalar(1.0));
  const Scalar _tmp117 = _tmp107 * _tmp108;
  const Scalar _tmp118 = _tmp117 * _tmp96 * (_tmp115 * _tmp48 + _tmp116 * _tmp65);
  const Scalar _tmp119 = -_tmp109 * _tmp18 + _tmp117 * _tmp28;
  const Scalar _tmp120 = Scalar(1.0) * _tmp101;
  const Scalar _tmp121 = _tmp65 * _tmp80;
  const Scalar _tmp122 = _tmp119 * _tmp96 * (-_tmp120 * _tmp121 + _tmp120 * _tmp48);
  const Scalar _tmp123 = Scalar(43.164000000000001) - fv1;
  const Scalar _tmp124 = _tmp84 * _tmp88 + _tmp87;
  const Scalar _tmp125 = 0;
  const Scalar _tmp126 = -_tmp110 * fh1 - _tmp118 * fh1 - _tmp122 * fh1 -
                         _tmp123 * _tmp96 * (-_tmp121 * _tmp125 + _tmp125 * _tmp48);
  const Scalar _tmp127 = Scalar(1.0) / (_tmp126);
  const Scalar _tmp128 = _tmp44 + _tmp85;
  const Scalar _tmp129 = _tmp128 * _tmp88;
  const Scalar _tmp130 = Scalar(1.0) / (-_tmp129 - _tmp40 + _tmp87);
  const Scalar _tmp131 = _tmp128 * _tmp130;
  const Scalar _tmp132 = _tmp103 * _tmp131 - _tmp104 * _tmp83 + _tmp92;
  const Scalar _tmp133 = Scalar(1.0) * _tmp86;
  const Scalar _tmp134 = Scalar(1.0) * _tmp130;
  const Scalar _tmp135 = Scalar(1.0) * _tmp109 * (_tmp103 * _tmp134 - _tmp132 * _tmp133);
  const Scalar _tmp136 = _tmp124 * _tmp130;
  const Scalar _tmp137 = _tmp86 * (-_tmp125 * _tmp83 - _tmp128 * _tmp136 + _tmp85);
  const Scalar _tmp138 = _tmp102 * _tmp134;
  const Scalar _tmp139 = _tmp128 * _tmp134;
  const Scalar _tmp140 = _tmp86 * (_tmp102 * _tmp139 - _tmp120 * _tmp83);
  const Scalar _tmp141 = Scalar(1.0) * _tmp119 * (_tmp138 - Scalar(1.0) * _tmp140);
  const Scalar _tmp142 = -_tmp112 + _tmp114 * _tmp131 - _tmp115 * _tmp83;
  const Scalar _tmp143 = Scalar(1.0) * _tmp117 * (_tmp114 * _tmp134 - _tmp133 * _tmp142);
  const Scalar _tmp144 = _tmp49 + _tmp53;
  const Scalar _tmp145 = _tmp144 * fh1;
  const Scalar _tmp146 = _tmp109 * _tmp145 + Scalar(5.1796800000000003) * _tmp24 + _tmp28 * fv1;
  const Scalar _tmp147 = _tmp139 * _tmp86;
  const Scalar _tmp148 = -Scalar(1.0) * _tmp134 + Scalar(1.0) * _tmp147;
  const Scalar _tmp149 = -_tmp117 * _tmp145 - Scalar(5.1796800000000003) * _tmp12 - _tmp18 * fv1;
  const Scalar _tmp150 = _tmp129 * _tmp134 + Scalar(1.0);
  const Scalar _tmp151 = _tmp134 * _tmp88;
  const Scalar _tmp152 = -Scalar(1.0) * _tmp133 * _tmp150 + Scalar(1.0) * _tmp151;
  const Scalar _tmp153 =
      Scalar(1.0) * _tmp123 * (-_tmp124 * _tmp134 - Scalar(1.0) * _tmp137 + Scalar(1.0)) +
      _tmp135 * fh1 + _tmp141 * fh1 + _tmp143 * fh1 + _tmp146 * _tmp148 + _tmp149 * _tmp152;
  const Scalar _tmp154 = std::asinh(_tmp127 * _tmp153);
  const Scalar _tmp155 = Scalar(1.0) * _tmp154;
  const Scalar _tmp156 = std::pow(_tmp126, Scalar(-2));
  const Scalar _tmp157 = -_tmp110 - _tmp118 - _tmp122;
  const Scalar _tmp158 = _tmp156 * _tmp157;
  const Scalar _tmp159 = _tmp109 * _tmp144;
  const Scalar _tmp160 = _tmp117 * _tmp144;
  const Scalar _tmp161 =
      (_tmp127 * (_tmp135 + _tmp141 + _tmp143 + _tmp148 * _tmp159 - _tmp152 * _tmp160) -
       _tmp153 * _tmp158) /
      std::sqrt(Scalar(std::pow(_tmp153, Scalar(2)) * _tmp156 + 1));
  const Scalar _tmp162 = Scalar(9.6622558468725703) * _tmp154;
  const Scalar _tmp163 =
      -_tmp126 * _tmp162 -
      Scalar(4.7744369927459998) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.2094487793051498) * _tmp70), Scalar(2)) +
                     Scalar(0.32387954179207445) *
                         std::pow(Scalar(1 - Scalar(0.36803241838814449) * _tmp73), Scalar(2))));
  const Scalar _tmp164 = Scalar(9.6622558468725703) * _tmp126;
  const Scalar _tmp165 = Scalar(0.1034955) * _tmp127;
  const Scalar _tmp166 = _tmp163 * _tmp165;
  const Scalar _tmp167 = _tmp109 * _tmp132 * _tmp86;
  const Scalar _tmp168 = _tmp119 * _tmp140;
  const Scalar _tmp169 = _tmp150 * _tmp86;
  const Scalar _tmp170 = _tmp117 * _tmp142 * _tmp86;
  const Scalar _tmp171 = _tmp134 * _tmp146;
  const Scalar _tmp172 = _tmp123 * _tmp137 - _tmp128 * _tmp171 * _tmp86 + _tmp149 * _tmp169 +
                         _tmp167 * fh1 + _tmp168 * fh1 + _tmp170 * fh1;
  const Scalar _tmp173 = _tmp105 * _tmp109;
  const Scalar _tmp174 = _tmp116 * _tmp117;
  const Scalar _tmp175 = _tmp119 * _tmp120;
  const Scalar _tmp176 = _tmp175 * fh1;
  const Scalar _tmp177 = _tmp123 * _tmp125;
  const Scalar _tmp178 = _tmp173 * fh1 + _tmp174 * fh1 - _tmp176 * _tmp80 - _tmp177 * _tmp80;
  const Scalar _tmp179 = Scalar(1.0) / (_tmp178);
  const Scalar _tmp180 = std::asinh(_tmp172 * _tmp179);
  const Scalar _tmp181 = Scalar(9.6622558468725703) * _tmp178;
  const Scalar _tmp182 =
      -_tmp180 * _tmp181 -
      Scalar(4.83288938413423) *
          std::sqrt(
              Scalar(std::pow(Scalar(1 - Scalar(0.20691555724053495) * _tmp62), Scalar(2)) +
                     Scalar(0.13818785160942856) *
                         std::pow(Scalar(-Scalar(0.55661923802822921) * _tmp59 - 1), Scalar(2))));
  const Scalar _tmp183 = Scalar(0.1034955) * _tmp179;
  const Scalar _tmp184 = _tmp182 * _tmp183;
  const Scalar _tmp185 = Scalar(1.0) * _tmp180;
  const Scalar _tmp186 = _tmp173 + _tmp174 - _tmp175 * _tmp80;
  const Scalar _tmp187 = Scalar(9.6622558468725703) * _tmp186;
  const Scalar _tmp188 = std::pow(_tmp178, Scalar(-2));
  const Scalar _tmp189 = _tmp186 * _tmp188;
  const Scalar _tmp190 = (-_tmp172 * _tmp189 + _tmp179 * (-_tmp147 * _tmp159 - _tmp160 * _tmp169 +
                                                          _tmp167 + _tmp168 + _tmp170)) /
                         std::sqrt(Scalar(std::pow(_tmp172, Scalar(2)) * _tmp188 + 1));
  const Scalar _tmp191 = _tmp104 * _tmp109;
  const Scalar _tmp192 = _tmp115 * _tmp117;
  const Scalar _tmp193 = _tmp176 + _tmp177 + _tmp191 * fh1 + _tmp192 * fh1;
  const Scalar _tmp194 = Scalar(1.0) / (_tmp193);
  const Scalar _tmp195 = _tmp119 * _tmp138;
  const Scalar _tmp196 = _tmp103 * _tmp109 * _tmp130;
  const Scalar _tmp197 = _tmp114 * _tmp117 * _tmp130;
  const Scalar _tmp198 = _tmp123 * _tmp136 - _tmp149 * _tmp151 + _tmp171 - _tmp195 * fh1 -
                         _tmp196 * fh1 - _tmp197 * fh1;
  const Scalar _tmp199 = std::asinh(_tmp194 * _tmp198);
  const Scalar _tmp200 = Scalar(1.0) * _tmp199;
  const Scalar _tmp201 = Scalar(9.6622558468725703) * _tmp193;
  const Scalar _tmp202 =
      -_tmp199 * _tmp201 -
      Scalar(8.3196563720703107) *
          std::sqrt(
              Scalar(Scalar(0.057067943527034905) *
                         std::pow(Scalar(-Scalar(0.50315118477277876) * _tmp41 - 1), Scalar(2)) +
                     std::pow(Scalar(-Scalar(0.12019727201198747) * _tmp45 - 1), Scalar(2))));
  const Scalar _tmp203 = Scalar(0.1034955) * _tmp194;
  const Scalar _tmp204 = _tmp202 * _tmp203;
  const Scalar _tmp205 = _tmp175 + _tmp191 + _tmp192;
  const Scalar _tmp206 = Scalar(9.6622558468725703) * _tmp205;
  const Scalar _tmp207 = std::pow(_tmp193, Scalar(-2));
  const Scalar _tmp208 = _tmp205 * _tmp207;
  const Scalar _tmp209 =
      (_tmp194 * (_tmp134 * _tmp159 + _tmp151 * _tmp160 - _tmp195 - _tmp196 - _tmp197) -
       _tmp198 * _tmp208) /
      std::sqrt(Scalar(std::pow(_tmp198, Scalar(2)) * _tmp207 + 1));

  // Output terms (1)
  Eigen::Matrix<Scalar, 4, 1> _res;

  _res(0, 0) = -Scalar(8.4690207536792048) * _tmp0 -
               Scalar(9.6622558468725703) * fh1 *
                   (-Scalar(1.0) * _tmp36 * _tmp37 * fv1 * std::sinh(_tmp38) -
                    Scalar(0.876505537412406) * _tmp36 -
                    (-Scalar(0.1034955) * _tmp33 * _tmp36 +
                     _tmp34 * (Scalar(9.6622558468725703) * _tmp30 * _tmp37 - _tmp32)) *
                        std::sinh(_tmp35)) +
               Scalar(9.6622558468725703) * std::cosh(_tmp35) -
               Scalar(9.6622558468725703) * std::cosh(_tmp38);
  _res(1, 0) =
      -Scalar(9.6622558468725703) * _tmp157 *
          (Scalar(0.86564762886483004) * _tmp127 + std::cosh(_tmp155) - std::cosh(_tmp166)) -
      _tmp164 *
          (-Scalar(0.86564762886483004) * _tmp158 + Scalar(1.0) * _tmp161 * std::sinh(_tmp155) -
           (-Scalar(0.1034955) * _tmp158 * _tmp163 +
            _tmp165 * (-_tmp157 * _tmp162 - _tmp161 * _tmp164)) *
               std::sinh(_tmp166));
  _res(2, 0) =
      -_tmp181 *
          (-Scalar(0.86627065637365697) * _tmp189 + Scalar(1.0) * _tmp190 * std::sinh(_tmp185) -
           (-Scalar(0.1034955) * _tmp182 * _tmp189 +
            _tmp183 * (-_tmp180 * _tmp187 - _tmp181 * _tmp190)) *
               std::sinh(_tmp184)) -
      _tmp187 * (Scalar(0.86627065637365697) * _tmp179 - std::cosh(_tmp184) + std::cosh(_tmp185));
  _res(3, 0) =
      -_tmp201 *
          (-Scalar(0.87679799777269396) * _tmp208 + Scalar(1.0) * _tmp209 * std::sinh(_tmp200) -
           (-Scalar(0.1034955) * _tmp202 * _tmp208 +
            _tmp203 * (-_tmp199 * _tmp206 - _tmp201 * _tmp209)) *
               std::sinh(_tmp204)) -
      _tmp206 * (Scalar(0.87679799777269396) * _tmp194 + std::cosh(_tmp200) - std::cosh(_tmp204));

  return _res;
}  // NOLINT(readability/fn_size)

// NOLINTNEXTLINE(readability/fn_size)
}  // namespace sym
