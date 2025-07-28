#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include"Levenshtein_distance.h"
#include "ecovlp.h"


cc_v scc = {0, 0, NULL, NULL};
cc_v scb = {0, 0, NULL, NULL};
cc_v sca = {0, 0, NULL, NULL};

// Function to create test data for cc_v structure
void create_test_cc_v(cc_v* x, size_t n_elements) {
    x->n = n_elements;
    x->m = n_elements;
    
    if (n_elements > 0) {
        x->a = (asg16_v*)calloc(n_elements, sizeof(asg16_v));
        x->f = (uint8_t*)calloc(n_elements, sizeof(uint8_t));
        
        for (size_t i = 0; i < n_elements; i++) {
            // Create test data for each asg16_v
            x->a[i].n = (i + 1) * 3; // Variable size
            x->a[i].m = x->a[i].n;
            
            if (x->a[i].n > 0) {
                x->a[i].a = (uint16_t*)calloc(x->a[i].n, sizeof(uint16_t));
                for (size_t j = 0; j < x->a[i].n; j++) {
                    x->a[i].a[j] = (uint16_t)(i * 100 + j); // Predictable test pattern
                }
            } else {
                x->a[i].a = NULL;
            }
            
            x->f[i] = (uint8_t)(i % 256); // Test flag values
        }
    } else {
        x->a = NULL;
        x->f = NULL;
    }
}

// Function to compare two cc_v structures
int compare_cc_v(const cc_v* a, const cc_v* b) {
    if (a->n != b->n || a->m != b->m) {
        printf("Size mismatch: a->n=%zu, b->n=%zu, a->m=%zu, b->m=%zu\n", 
               a->n, b->n, a->m, b->m);
        return 0;
    }
    
    if (a->n == 0) {
        return (a->a == NULL && b->a == NULL && a->f == NULL && b->f == NULL);
    }
    
    if ((a->a == NULL) != (b->a == NULL) || (a->f == NULL) != (b->f == NULL)) {
        printf("Null pointer mismatch\n");
        return 0;
    }
    
    // Compare asg16_v arrays
    for (size_t i = 0; i < a->n; i++) {
        if (a->a[i].n != b->a[i].n || a->a[i].m != b->a[i].m) {
            printf("asg16_v[%zu] size mismatch: a->n=%zu, b->n=%zu, a->m=%zu, b->m=%zu\n", 
                   i, a->a[i].n, b->a[i].n, a->a[i].m, b->a[i].m);
            return 0;
        }
        
        if (a->a[i].n > 0) {
            if ((a->a[i].a == NULL) != (b->a[i].a == NULL)) {
                printf("asg16_v[%zu] null pointer mismatch\n", i);
                return 0;
            }
            
            if (a->a[i].a != NULL && memcmp(a->a[i].a, b->a[i].a, a->a[i].n * sizeof(uint16_t)) != 0) {
                printf("asg16_v[%zu] data mismatch\n", i);
                return 0;
            }
        }
    }
    
    // Compare flag arrays
    if (memcmp(a->f, b->f, a->n * sizeof(uint8_t)) != 0) {
        printf("Flag array mismatch\n");
        return 0;
    }
    
    return 1;
}

// Function to cleanup cc_v structure
void cleanup_test_cc_v(cc_v* x) {
    if (x->a != NULL) {
        for (size_t i = 0; i < x->n; i++) {
            if (x->a[i].a != NULL) {
                free(x->a[i].a);
            }
        }
        free(x->a);
        x->a = NULL;
    }
    if (x->f != NULL) {
        free(x->f);
        x->f = NULL;
    }
    x->n = x->m = 0;
}

// Test function for individual cc_v read/write
int test_cc_v_io() {
    printf("Testing individual cc_v read/write...\n");
    
    cc_v original = {0, 0, NULL, NULL};
    cc_v loaded = {0, 0, NULL, NULL};
    
    // Test with various sizes including edge cases
    size_t test_sizes[] = {0, 1, 5, 100};
    size_t num_tests = sizeof(test_sizes) / sizeof(test_sizes[0]);
    
    for (size_t t = 0; t < num_tests; t++) {
        printf("  Testing with %zu elements...\n", test_sizes[t]);
        
        // Create test data
        create_test_cc_v(&original, test_sizes[t]);
        
        // Write to file
        FILE* fp = fopen("test_cc_v.bin", "wb");
        if (fp == NULL) {
            printf("ERROR: Cannot create test file\n");
            cleanup_test_cc_v(&original);
            return 0;
        }
        write_cc_v(&original, fp);
        fclose(fp);
        
        // Load from file
        fp = fopen("test_cc_v.bin", "rb");
        if (fp == NULL) {
            printf("ERROR: Cannot open test file for reading\n");
            cleanup_test_cc_v(&original);
            return 0;
        }
        if (!load_cc_v(&loaded, fp)) {
            printf("ERROR: Failed to load cc_v data\n");
            fclose(fp);
            cleanup_test_cc_v(&original);
            return 0;
        }
        fclose(fp);
        
        // Compare
        if (!compare_cc_v(&original, &loaded)) {
            printf("ERROR: Data mismatch for size %zu\n", test_sizes[t]);
            cleanup_test_cc_v(&original);
            cleanup_test_cc_v(&loaded);
            unlink("test_cc_v.bin");
            return 0;
        }
        
        printf("  ✓ Size %zu test passed\n", test_sizes[t]);
        
        // Cleanup
        cleanup_test_cc_v(&original);
        cleanup_test_cc_v(&loaded);
        unlink("test_cc_v.bin");
    }
    
    return 1;
}

// Test function for the full cc_v_all read/write
int test_cc_v_all_io() {
    printf("Testing cc_v_all read/write...\n");

    
    // Create test data for global structures
    create_test_cc_v(&scc, 10);
    create_test_cc_v(&scb, 15);
    create_test_cc_v(&sca, 8);
    
    // Write all structures
    write_cc_v_all("test_output");
    
    // Backup the test data
    cc_v scc_test = scc;
    cc_v scb_test = scb;
    cc_v sca_test = sca;
    
    // Clear global structures
    memset(&scc, 0, sizeof(cc_v));
    memset(&scb, 0, sizeof(cc_v));
    memset(&sca, 0, sizeof(cc_v));
    
    // Load all structures
    if (!load_cc_v_all("test_output")) {
        printf("ERROR: Failed to load cc_v_all data\n");
        // Restore original structures
       
        return 0;
    }
    
    // Compare loaded data with original
    int success = 1;
    if (!compare_cc_v(&scc, &scc_test)) {
        printf("ERROR: scc data mismatch\n");
        success = 0;
    }
    if (!compare_cc_v(&scb, &scb_test)) {
        printf("ERROR: scb data mismatch\n");
        success = 0;
    }
    if (!compare_cc_v(&sca, &sca_test)) {
        printf("ERROR: sca data mismatch\n");
        success = 0;
    }
    
    if (success) {
        printf("  ✓ cc_v_all test passed\n");
    }
    
    // Cleanup test data
    cleanup_test_cc_v(&scc);
    cleanup_test_cc_v(&scb);
    cleanup_test_cc_v(&sca);
    cleanup_test_cc_v(&scc_test);
    cleanup_test_cc_v(&scb_test);
    cleanup_test_cc_v(&sca_test);
    
    // Restore original structures
   
    
    // Cleanup test files
    unlink("test_output.scc.bin");
    unlink("test_output.scb.bin");
    unlink("test_output.sca.bin");
    
    return success;
}

int main() {
    printf("Running cc_v I/O unit tests...\n\n");
    
    int all_passed = 1;
    
    // Test individual cc_v I/O
    if (!test_cc_v_io()) {
        printf("❌ Individual cc_v I/O test FAILED\n");
        all_passed = 0;
    } else {
        printf("✅ Individual cc_v I/O test PASSED\n");
    }
    
    printf("\n");
    
    // Test cc_v_all I/O
    if (!test_cc_v_all_io()) {
        printf("❌ cc_v_all I/O test FAILED\n");
        all_passed = 0;
    } else {
        printf("✅ cc_v_all I/O test PASSED\n");
    }
    
    printf("\n");
    
    if (all_passed) {
        printf("🎉 All tests PASSED! cc_v I/O functions preserve data integrity.\n");
        return 0;
    } else {
        printf("💥 Some tests FAILED! cc_v I/O functions have issues.\n");
        return 1;
    }
}
