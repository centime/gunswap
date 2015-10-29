module.exports = function(grunt)
{
    // Configuration.
    grunt.initConfig({

        // Import package information.
        package: grunt.file.readJSON('package.json'),

        browserify: 
        {
            main: {
                src: 'js/index.js',
                dest: 'js/build/gunswap.js'
            }            
        },

        // Test.
        simplemocha:
        {
            options: {
                timeout: 3000,
                ignoreLeaks: false,
                ui: 'bdd',
                reporter: 'nyan'
            },

            all:
            {
                src: [ 'test.js' ]
            }
        },

        uglify: 
        {
            dist: {
                files: {
                    'js/build/gunswap.min.js': ['js/build/gunswap.js']
                }
            }
        },

        // Watch.
        watch:
        {
            // Build product.
            default:
            {
                files: ['**'],
                tasks: ['build'],
            }

        }

    });

    // Load plugins.
    grunt.loadNpmTasks('grunt-simple-mocha');
    grunt.loadNpmTasks('grunt-contrib-watch');
    grunt.loadNpmTasks('grunt-contrib-concat');
    grunt.loadNpmTasks('grunt-contrib-uglify');
    grunt.loadNpmTasks('grunt-browserify');

    // `default` task for the everyday.
    grunt.registerTask('default', ['watch:default']);

    // `build`, `test` task for Travis CI.
    grunt.registerTask('test', ['simplemocha']);

    grunt.registerTask('build', ['browserify','uglify']);

};