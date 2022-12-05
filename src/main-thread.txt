#include <iostream>
#include <vector>
#include <functional> // for std::reference_wrapper<T>
#include <future> // for std::async
#include <mutex> // for std::mutex
#include <chrono> // for Timer
#include <thread> // for std::this_thread::sleep_for(seconds(5));

using namespace std::chrono;

struct Mesh;
class LoadMesh;
struct Timer;

// Timer class -> used for tracking the run time
struct Timer //for counting the time
{
    std::chrono::time_point<std::chrono::steady_clock>start, end;
    std::chrono::duration<float>duration;

    Timer() //set default value
    {
        start = end = std::chrono::high_resolution_clock::now();
        duration = end - start;
    }

    ~Timer() // get the end value and print the duration time
    {
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;

        std::cout << "Time: " << duration.count() << "s\n";
    }
};


template<typename T>
using Ref = std::reference_wrapper<T>;


struct Mesh {
    static Ref<Mesh> Load(const std::string& file) {
        // ...
        std::this_thread::sleep_for(seconds(5));
        Mesh m;
        return std::ref(m); // can not use std::ref(Mesh())

        // generate an object of type std::reference_wrapper
    }
};


namespace GenerateMesh {

    std::vector<Ref<Mesh>> m_Meshes;
    std::vector<std::future<void>> m_Futures; // store the return value of std::async, necessary step to make async work
    std::mutex s_MeshesMutex;


    void LoadMesh(std::vector<Ref<Mesh>>* meshes, std::string filepath /* make copy */) {
        auto mesh = Mesh::Load(filepath);

        // using a local lock_guard to lock mtx guarantees unlocking on destruction / exception:
        std::lock_guard<std::mutex> lock(s_MeshesMutex); // lock the meshes to avoid conflict
        meshes->emplace_back(mesh);

        //std::cout<<mesh;
    }


    void LoadMeshes() {
        std::vector<std::string> meshFilepaths; // read from a .txt file ...

        // ...
        // meshFilepaths.emplace_back(...)
        for (int i = 0; i != 10; ++i) {
            meshFilepaths.emplace_back("");
        }

        std::cout << "size of meshFilepaths: " << meshFilepaths.size() << std::endl;

        /*for(const auto& file : meshFilepaths){
            m_Meshes.emplace_back(Mesh::Load(file));
        }*/

        /*
        * it is important to save the result of std::async()
        * to enable the async process
        */
        for (const auto& file : meshFilepaths) {
            m_Futures.emplace_back(std::async(std::launch::async, LoadMesh, &m_Meshes, file));
        }


        /*
        * if we wish to get the result value and keep processing
        * we need to use get() of every future object
        */
        for (auto& futureObject : m_Futures) {
            futureObject.get();
        }
    }

};


int main()
{
    Timer timer;
    GenerateMesh::LoadMeshes();
    std::cout << "done" << '\n';
    return 0;
}