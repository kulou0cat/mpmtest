#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>

constexpr int dim = 2;//二维模拟
constexpr int n_particles = 8192;//粒子数量
constexpr int n_grid = 128;//网格数量
constexpr float dx = 1.0f / n_grid;// 网格间距
constexpr float inv_dx = 1.0f / dx;// 1/dx
constexpr float dt = 2.0e-4f;// 时间步长
constexpr float p_vol = (dx * 0.5f) * (dx * 0.5f);// 粒子体积
constexpr float p_rho = 1.0f;// 粒子密度
constexpr float p_mass = p_vol * p_rho;// 粒子质量
constexpr float E = 400.0f;// 弹性模量

// 粒子属性
std::vector<Eigen::Vector2f> x(n_particles);// 粒子位置
std::vector<Eigen::Vector2f> v(n_particles);// 粒子速度
std::vector<Eigen::Matrix2f> C(n_particles, Eigen::Matrix2f::Zero());// 粒子仿射矩阵Cn
std::vector<float> J(n_particles, 1.0f);// 粒子体积变化率

// 网格属性 (按行展平存储: index = i * n_grid + j)
std::vector<Eigen::Vector2f> grid_v(n_grid* n_grid, Eigen::Vector2f::Zero());// 网格速度
std::vector<float> grid_m(n_grid* n_grid, 0.0f);// 网格质量
inline int grid_idx(int i, int j) {
    return i * n_grid + j;
}

void substep() {
    // 1. 重置网格数据,全部置为0
    std::fill(grid_v.begin(), grid_v.end(), Eigen::Vector2f::Zero());
    std::fill(grid_m.begin(), grid_m.end(), 0.0f);

    // 2. P2G (Particle to Grid: 粒子到网格传输)
    for (int p = 0; p < n_particles; ++p) {
        Eigen::Vector2i base = (x[p] * inv_dx - Eigen::Vector2f::Constant(0.5f)).cast<int>();
        Eigen::Vector2f fx = x[p] * inv_dx - base.cast<float>();

        // 二次 B 样条权重 (Quadratic B-spline weights)
        Eigen::Vector2f w[3] = {
            0.5f * (Eigen::Vector2f::Constant(1.5f) - fx).array().square().matrix(),//w0
            Eigen::Vector2f::Constant(0.75f) - (fx - Eigen::Vector2f::Constant(1.0f)).array().square().matrix(),//w1
            0.5f * (fx - Eigen::Vector2f::Constant(0.5f)).array().square().matrix()//w2
        };

		float stress = -dt * p_vol * (J[p] - 1.0f) * 4.0f * inv_dx * inv_dx * E;// 弹性应力
		Eigen::Matrix2f affine = stress * Eigen::Matrix2f::Identity() + p_mass * C[p];// 仿射矩阵 Ap=SpI+mpCp=内部应力+APIC仿射速度

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                int gx = base.x() + i;
                int gy = base.y() + j;
				if (gx < 0 || gx >= n_grid || gy < 0 || gy >= n_grid) continue;// 越界检查

				Eigen::Vector2f dpos = (Eigen::Vector2f(static_cast<float>(i), static_cast<float>(j)) - fx) * dx;// 粒子到网格节点的偏移量xi-xp
                float weight = w[i].x() * w[j].y();

				int idx = grid_idx(gx, gy);// 网格节点索引
				grid_v[idx] += weight * (p_mass * v[p] + affine * dpos);// 网格动量更新Ev
				grid_m[idx] += weight * p_mass;// 网格质量更新
            }
        }
    }

    // 3. 网格动量更新与边界条件处理
    constexpr int bound = 3;

    for (int i = 0; i < n_grid; ++i) {
        for (int j = 0; j < n_grid; ++j) {
            int idx = grid_idx(i, j);
            if (grid_m[idx] > 0.0f) {
				float inv_m = 1.0f / grid_m[idx];// 网格质量的倒数
                //计算网格速度
				grid_v[idx] *= inv_m;// 网格速度更新 Ev=mv  v=Ev/m
				grid_v[idx].y() -= dt * 9.8f;// 重力加速度
				// 边界条件处理，防止粒子穿过边界
                if (i < bound && grid_v[idx].x() < 0.0f) grid_v[idx].x() = 0.0f;
                if (i > n_grid - bound && grid_v[idx].x() > 0.0f) grid_v[idx].x() = 0.0f;
                if (j < bound && grid_v[idx].y() < 0.0f) grid_v[idx].y() = 0.0f;
                if (j > n_grid - bound && grid_v[idx].y() > 0.0f) grid_v[idx].y() = 0.0f;
            }
        }
    }

    // 4. G2P (Grid to Particle: 网格到粒子传输)
    for (int p = 0; p < n_particles; ++p) {
        Eigen::Vector2i base = (x[p] * inv_dx - Eigen::Vector2f::Constant(0.5f)).cast<int>();
        Eigen::Vector2f fx = x[p] * inv_dx - base.cast<float>();

        Eigen::Vector2f w[3] = {
            0.5f * (Eigen::Vector2f::Constant(1.5f) - fx).array().square().matrix(),
            Eigen::Vector2f::Constant(0.75f) - (fx - Eigen::Vector2f::Constant(1.0f)).array().square().matrix(),
            0.5f * (fx - Eigen::Vector2f::Constant(0.5f)).array().square().matrix()
        };

		Eigen::Vector2f new_v = Eigen::Vector2f::Zero();// 新的粒子速度
		Eigen::Matrix2f new_C = Eigen::Matrix2f::Zero();// 新的粒子仿射矩阵Cn

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                int gx = base.x() + i;//x
                int gy = base.y() + j;//y
				if (gx < 0 || gx >= n_grid || gy < 0 || gy >= n_grid) continue;// 越界检查

                Eigen::Vector2f dpos = Eigen::Vector2f(static_cast<float>(i), static_cast<float>(j)) - fx;// 粒子到网格节点的偏移量(xi-xp)/dx
				const Eigen::Vector2f& g_v = grid_v[grid_idx(gx, gy)];//获取网格速度
                float weight = w[i].x() * w[j].y();

                new_v += weight * g_v;
                // outer_product: g_v * dpos^T
                new_C += 4.0f * weight * (g_v * dpos.transpose()) * inv_dx;
            }
        }

        v[p] = new_v;
        x[p] += dt * v[p];
        J[p] *= 1.0f + dt * new_C.trace();
        C[p] = new_C;
    }
}


int main() {
    constexpr int win_size = 512;
    constexpr int substeps_per_frame = 50;
    const char* window_name = "MPM88";
    // ==========================================
    // 1. 鼠标交互状态
    // ==========================================
    struct MouseState {
        bool left_down = false;
        cv::Point pos{ 0, 0 };
        cv::Point previous_pos{ 0, 0 };
    };
    MouseState mouse;
    // ==========================================
    // 2. 初始化 / 重置模拟
    // ==========================================
    auto reset_simulation = [&]() {
        // 固定随机种子，保证每次 R 重置结果相同
        std::mt19937 rng(42);
        std::uniform_real_distribution<float> dist(0.0f, 1.0f);
        for (int p = 0; p < n_particles; ++p) {
            x[p] = Eigen::Vector2f(
                dist(rng) * 0.4f + 0.2f,
                dist(rng) * 0.4f + 0.2f
            );

            v[p] = Eigen::Vector2f(0.0f, -1.0f);
            C[p].setZero();
            J[p] = 1.0f;
        }
        std::fill(
            grid_v.begin(),
            grid_v.end(),
            Eigen::Vector2f::Zero()
        );
        std::fill(grid_m.begin(), grid_m.end(), 0.0f);
        mouse.left_down = false;
        mouse.previous_pos = mouse.pos;

        std::cout << "Simulation reset." << std::endl;
        };
    reset_simulation();
    // ==========================================
    // 3. 创建 OpenCV 窗口
    // ==========================================
    cv::namedWindow(window_name, cv::WINDOW_AUTOSIZE);
    // 注册鼠标事件回调
    cv::setMouseCallback(
        window_name,
        [](int event, int mx, int my, int flags, void* userdata) {
            auto* state = static_cast<MouseState*>(userdata);

            if (event == cv::EVENT_LBUTTONDOWN) {
                state->left_down = true;
                state->pos = cv::Point(mx, my);
                state->previous_pos = state->pos;
            }
            else if (event == cv::EVENT_MOUSEMOVE) {
                state->pos = cv::Point(mx, my);
            }
            else if (event == cv::EVENT_LBUTTONUP) {
                state->left_down = false;
                state->pos = cv::Point(mx, my);
                state->previous_pos = state->pos;
            }
        },
        &mouse
    );
    // ==========================================
    // 4. 创建画布
    // ==========================================
    cv::Mat canvas(win_size, win_size, CV_8UC3);
    const cv::Scalar bg_color(65, 47, 17);
    const cv::Scalar particle_color(135, 133, 6);
    std::cout << "MPM88 Simulation started." << std::endl;
    std::cout << "R: Reset simulation" << std::endl;
    std::cout << "Left Mouse: Push / Stir particles" << std::endl;
    std::cout << "ESC: Exit" << std::endl;
    // ==========================================
    // 5. 主模拟循环
    // =========================================
    for (int frame = 0; frame < 20000; ++frame) {

        // ======================================
        // 5.1 鼠标与粒子交互
        // ======================================

        // ======================================
        // 5.1 鼠标与粒子交互（增强版本）
        // ======================================

        if (mouse.left_down) {

            // 鼠标位置转换到 MPM 坐标
            Eigen::Vector2f mouse_pos(
                static_cast<float>(mouse.pos.x) / win_size,
                1.0f - static_cast<float>(mouse.pos.y) / win_size
            );

            // 鼠标拖动方向
            Eigen::Vector2f mouse_drag(
                static_cast<float>(
                    mouse.pos.x - mouse.previous_pos.x
                    ) / win_size,

                -static_cast<float>(
                    mouse.pos.y - mouse.previous_pos.y
                    ) / win_size
            );

            // 限制鼠标单帧移动距离
            constexpr float max_drag = 0.08f;

            float drag_length = mouse_drag.norm();

            if (drag_length > max_drag) {
                mouse_drag *= max_drag / drag_length;
            }

            // ==================================
            // 鼠标交互参数
            // ==================================

            // 影响半径：约 72 像素
            constexpr float mouse_radius = 0.14f;

            // 点击推动强度
            constexpr float push_strength = 2.5f;

            // 鼠标拖动强度
            constexpr float drag_strength = 60.0f;

            // 限制单帧速度冲量
            constexpr float max_impulse = 5.0f;

            // ==================================
            // 对附近粒子施加速度冲量
            // ==================================

            for (int p = 0; p < n_particles; ++p) {

                Eigen::Vector2f offset = x[p] - mouse_pos;

                float distance = offset.norm();

                if (distance >= mouse_radius) {
                    continue;
                }

                // 线性距离衰减
                float influence =
                    1.0f - distance / mouse_radius;

                // 径向推动方向
                Eigen::Vector2f push_direction =
                    Eigen::Vector2f::Zero();

                if (distance > 1.0e-6f) {
                    push_direction = offset / distance;
                }

                // 点击产生径向排斥
                Eigen::Vector2f push_velocity =
                    push_strength * push_direction;

                // 拖动产生定向速度
                Eigen::Vector2f drag_velocity =
                    drag_strength * mouse_drag;

                // 合成鼠标作用
                Eigen::Vector2f impulse =
                    influence * (push_velocity + drag_velocity);

                // 限制冲量，避免数值不稳定
                float impulse_length = impulse.norm();

                if (impulse_length > max_impulse) {
                    impulse *= max_impulse / impulse_length;
                }

                // 将冲量施加到粒子速度
                v[p] += impulse;
            }
        }

        // 更新鼠标位置
        mouse.previous_pos = mouse.pos;

        // ======================================
        // 5.2 MPM 物理模拟
        // ======================================
        for (int s = 0; s < substeps_per_frame; ++s) {
            substep();
        }
        // ======================================
        // 5.3 渲染
        // =====================================
        canvas.setTo(bg_color);
        for (int p = 0; p < n_particles; ++p) {

            int px = static_cast<int>(
                x[p].x() * win_size
                );
            int py = static_cast<int>(
                (1.0f - x[p].y()) * win_size
                );
            if (px >= 0 && px < win_size &&
                py >= 0 && py < win_size) {
                cv::circle(
                    canvas,
                    cv::Point(px, py),
                    2,
                    particle_color,
                    cv::FILLED,
                    cv::LINE_AA
                );
            }
        }
        cv::imshow(window_name, canvas);

        // ======================================
        // 5.4 键盘事件
        // ======================================
        int key = cv::waitKey(1);
        // ESC：退出
        if (key == 27) {
            break;
        }
        // R：重置模拟
        if (key == 'r' || key == 'R') {
            reset_simulation();

            // 重新开始帧计数
            frame = -1;
            continue;
        }
    }

    cv::destroyAllWindows();
    return 0;
}
