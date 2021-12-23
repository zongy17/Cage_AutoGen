#include "glRender.h"


static void framebuffer_size_callback(GLFWwindow* window, int width, int height);
static void mouse_callback(GLFWwindow* window, double xpos, double ypos);
static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
static void processInput(GLFWwindow *window);

// settings
static const unsigned int SCR_WIDTH = 800;
static const unsigned int SCR_HEIGHT = 600;

// camera
static Camera camera(glm::vec3(2.0f, 2.0f, 1.0f));

// 鼠标状态
static bool CURSOR_DISABLED = true;//是否选择cursor追随模式（FPS相机）
static bool MOUSE_LEFT_ON = false;
static bool MOUSE_RIGHT_ON = false;
static float lastX = SCR_WIDTH / 2.0f;
static float lastY = SCR_HEIGHT / 2.0f;
static bool firstMouse = true;

// timing
static float deltaTime = 0.0f;
static float lastFrame = 0.0f;

// shader files
static std::string mesh_vs_file = "../resource/mesh.vs";
static std::string mesh_fs_file = "../resource/mesh.fs";
static std::string cage_vs_file = "../resource/cage.vs";
static std::string cage_fs_file = "../resource/cage.fs";

// window
static GLFWwindow* window;
int cagePtr_id = 1;
zyMesh * cagePtr = nullptr;

void renderLoop(std::string mesh_file, std::string voxel_file, std::string cage_file) {
    // glfw: initialize and configure
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);//如果需要请求调试输出

    // glfw window creation
    window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Cage automatic generation", NULL, NULL);
    if (window == NULL){
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        exit(1);
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    if (CURSOR_DISABLED) {
        glfwSetCursorPosCallback(window, mouse_callback);
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);// tell GLFW to capture our mouse
    }
    glfwSetScrollCallback(window, scroll_callback);

    // glad: load all OpenGL function pointers
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)){
        std::cout << "Failed to initialize GLAD" << std::endl;
        exit(1);
    }

    // configure global opengl state
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    // build and compile shaders
    Shader meshShader(mesh_vs_file.c_str(), mesh_fs_file.c_str());
    Shader cageShader(cage_vs_file.c_str(), cage_fs_file.c_str());
    // load meshes
    zyMesh ourMesh(mesh_file, false, true);
    zyMesh ourVoxl(voxel_file,false, true);
    zyMesh ourCage(cage_file, false, true);
    cagePtr = & ourCage;

    glLineWidth(1.0f);
    glPointSize(2.0f);

    // glEnable(GL_CULL_FACE);//面剔除，如果想要背面也画出来（尤其是GL_LINE和GL_POINT模式下），则关掉这两行
    // glCullFace(GL_BACK);//剔除掉背向面

    // render loop
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        processInput(window);
        if (cagePtr_id == 1)      cagePtr = & ourCage;
        else if (cagePtr_id == 2) cagePtr = & ourVoxl;

        // render
        glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // view/projection transformations
        glm::mat4 model;
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        glm::mat4 view = camera.GetViewMatrix();
        

        // 渲染原始的mesh网格
        meshShader.use();// don't forget to enable shader before setting uniforms
        meshShader.setMat4("projection", projection);
        meshShader.setMat4("view", view);
        meshShader.setMat4("model", glm::mat4(1.0f));
        // model = glm::mat4(1.0f);
        // model = glm::translate(model, glm::vec3(0.0f, 0.0f, 0.0f)); // translate it down so it's at the center of the scene
        // model = glm::scale(model, glm::vec3(1.0f, 1.0f, 1.0f));	// it's a bit too big for our scene, so scale it down
        // meshShader.setMat4("model", model);
        
        meshShader.setBool("sgColor", true);
        meshShader.setVec3("backColor", glm::vec3(0.8f, 0.0f, 0.0f));

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);//面
        // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);//边（网格线 wireframe）
        // glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);//点
        ourMesh.Draw(meshShader);

        // 面+边（点）模式，则先画面，再画边（点）
        meshShader.setBool("sgColor", true);
        meshShader.setVec3("backColor", glm::vec3(0.0f, 0.0f, 0.0f));
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);//在面的基础上叠加上边
        glEnable(GL_POLYGON_OFFSET_LINE);
        // glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);//在面的基础上叠加上点
        // glEnable(GL_POLYGON_OFFSET_POINT);
        glPolygonOffset(-0.5, -0.5); // move closer to camera
        ourMesh.Draw(meshShader);

        // 渲染光滑后的cage网格
        cageShader.use();
        cageShader.setMat4("projection", projection);
        cageShader.setMat4("view", view);
        cageShader.setMat4("model", glm::mat4(1.0f));

        cageShader.setBool("sgColor", true);
        cageShader.setVec3("backColor", glm::vec3(0.0f, 1.0f, 0.0f));

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);//面
        // ourCage.Draw(cageShader);
        cagePtr->Draw(cageShader);

        // 面+边（点）模式，则先画面，再画边（点）
        cageShader.setBool("sgColor", true);
        cageShader.setVec3("backColor", glm::vec3(0.0f, 0.0f, 0.0f));// 黑色
        glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);//在面的基础上叠加上边
        glEnable(GL_POLYGON_OFFSET_LINE);
        glPolygonOffset(-0.5, -0.5); // move closer to camera
        // ourCage.Draw(cageShader);
        cagePtr->Draw(cageShader);

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // glfw: terminate, clearing all previously allocated GLFW resources.
    glfwTerminate();
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);


    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {//正在被按下
        if (MOUSE_LEFT_ON == false) {//之前没被按下
            MOUSE_LEFT_ON = true;
            //保存下来lastX和lastY
            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);
            lastX = xpos;
            lastY = ypos;
        } else {//之前一直被按着
            //重新获得当前的位置
            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);
            float xoffset = xpos - lastX;
            float yoffset = lastY - ypos;
            lastX = xpos;
            lastY = ypos;

            camera.ProcessMouseMovement(-xoffset, -yoffset);//更改摄像机位置
        }
    } else {//被松开
        MOUSE_LEFT_ON = false;
    }

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {//正在被按下
        if (MOUSE_RIGHT_ON == false) {//之前没被按下
            MOUSE_RIGHT_ON = true;
            //保存下来lastX和lastY
            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);
            lastX = xpos;
            lastY = ypos;
        } else {//之前一直被按着
            //重新获得当前的位置
            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);
            float sensitivity = 3e-3;
            float xoffset = (xpos - lastX)*sensitivity;
            float yoffset = (lastY - ypos)*sensitivity;
            lastX = xpos;
            lastY = ypos;

            camera.ProcessKeyboard(RIGHT, -xoffset);
            camera.ProcessKeyboard(UP   , -yoffset);
        }
    } else {//被松开
        MOUSE_RIGHT_ON = false;
    }
    
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
        cagePtr_id = 1;
    else if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
        cagePtr_id = 2;

}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

// glfw: whenever the mouse moves, this callback is called
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);
}