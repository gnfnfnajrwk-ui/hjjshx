# PID Tuning Simulator

브라우저에서 바로 실행되는 PID 튜닝 시뮬레이터입니다. 서버가 필요 없으며 GitHub Pages로 자동 배포할 수 있습니다.

## 기본 공정

현재 기본 공정은 다음 FOPDT 모델입니다.

```text
300 · dT/dt = -(T - 25) + 16 · u(t - 30)
```

- 초기 온도: 25 ℃
- 설정값: 1000 ℃
- 시간지연: 30 s
- 출력 제한: 0~100 %
- Euler 적분, dt=1 s
- 총 시뮬레이션: 5400 s
- 기본 PID: Kc=0.3, Ti=200 s, Td=25 s
- 출력 포화 시 적분 누적 중지(anti-windup)

## 주요 기능

- Kc, Ti, Td 슬라이더/숫자 입력을 바꾸면 즉시 재계산
- FOPDT 파라미터를 UI에서 자유롭게 변경
- Custom ODE 모드에서 `dy/dt` 식을 직접 입력
- 설정값, 시뮬레이션 시간, dt, 출력 한계 변경
- 현재 PID와 기본 PID 응답 비교
- 최종값, 최종오차, 오버슈트, ±2% 정착시간, IAE, 포화시간 표시
- 공정값/제어출력 그래프
- CSV 내보내기
- 모든 계산은 브라우저 안에서 실행

## Custom ODE

지원 변수:

- `y`: 현재 공정값
- `u`: 시간지연이 적용된 제어 출력
- `t`: 현재 시간 [s]
- `sp`: 설정값

지원 연산자: `+ - * / ^`

지원 함수: `abs`, `sin`, `cos`, `tan`, `exp`, `log`, `sqrt`, `min`, `max`, `pow`

예시:

```text
(-(y - 25) + 16*u) / 300
```

## 로컬 실행

`index.html`을 브라우저에서 열면 됩니다. 더 안정적인 로컬 서버를 원하면:

```bash
python -m http.server 8080
```

그 후 `http://localhost:8080` 접속.

## GitHub Pages 자동 배포

`main` 브랜치에 push하면 `.github/workflows/pages.yml`이 실행되어 GitHub Pages로 배포됩니다.

Repository Settings → Pages에서 Source가 **GitHub Actions**로 설정되어 있어야 합니다.

## License

MIT
